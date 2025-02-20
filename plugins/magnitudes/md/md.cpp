/************************************************************************
 *                                                                      *
 * Copyright (C) 2012 OVSM/IPGP                                         *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                  *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * This program is part of 'Projet TSUAREG - INTERREG IV Caraïbes'.     *
 * It has been co-financed by the European Union and le Ministère de    *
 * l'Ecologie, du Développement Durable, des Transports et du Logement. *
 *                                                                      *
 ************************************************************************/


#define SEISCOMP_COMPONENT Md

#include "md.h"
#include "l4c1hz.h"
#include <seiscomp/processing/waveformprocessor.h>
#include <seiscomp/math/filter/stalta.h>
#include <seiscomp/math/filter/butterworth.h>
#include <seiscomp/logging/log.h>
#include <seiscomp/core/strings.h>
#include <seiscomp/seismology/magnitudes.h>
#include <seiscomp/math/mean.h>
#include <seiscomp/math/filter/seismometers.h>
#include <seiscomp/math/restitution/fft.h>
#include <seiscomp/math/geo.h>
#include <iostream>
#include <unistd.h>


#define _DEPTH_MAX 200.0
#define _SIGNAL_WINDOW_END 30.0
#define _LINEAR_CORRECTION 1.0
#define _OFFSET 0.0
#define _SNR_MIN 1.2
#define _TAPER 5
#define _SIGNAL_LENGTH 150
#define _DELTA_MAX 400.0
#define _MD_MAX 5.0
#define _FMA -0.87
#define _FMB 2.0
#define _FMD 0.0035
#define _FMF 0.0
#define _FMZ 0.0
#define _STACOR 0.0 //! not fully implemented !
/**
 * Seismometer selection
 * 1 for Wood-Anderson
 * 2 for Seismometer 5 sec
 * 3 for WWSSN LP ? filter
 * 4 for WWSSN SP? filter
 * 5 for Generic Seismometer ? filter
 * 6 for Butterworth Low Pass ? filter
 * 7 for Butterwoth High Pass ? filter
 * 8 for Butterworth Band Pass ? filter
 * 9 for L4C 1Hz seismometer
 **/
#define _SEISMO 9
#define _BUTTERWORTH ""

ADD_SC_PLUGIN("Md duration magnitude plugin", "IPGP <www.ipgp.fr>", 0, 1, 2)

#define AMPTAG "[Amp] [Md]"
#define MAGTAG "[Mag] [Md]"


using namespace Seiscomp;
using namespace Seiscomp::Math;
using namespace Seiscomp::Processing;


namespace {


template <typename T>
T getParam(const Settings &settings, const char *name) {
	T v;

	if ( settings.getValue(v, (std::string("magnitudes.") + name).c_str()) )
		return v;

	if ( settings.getValue(v, name) ) {
		SEISCOMP_WARNING("Configure magnitudes.%s in global bindings. "
		                 "The old parameter %s has been deprecated "
		                 "and should be replaced.", name, name);
		return v;
	}

	throw Seiscomp::Config::OptionNotFoundException(name);
}


}


/*----[ AMPLITUDE PROCESSOR CLASS ]----*/

REGISTER_AMPLITUDEPROCESSOR(AmplitudeProcessor_Md, "Md");

struct ampConfig {

		double DEPTH_MAX;
		double SIGNAL_WINDOW_END;
		double SNR_MIN;
		double TAPER;
		double SIGNAL_LENGTH;
		double DELTA_MAX;
		double MD_MAX;
		double FMA;
		double FMB;
		double FMD;
		double FMF;
		double FMZ;
		double STACOR;
		int SEISMO;
		std::string BUTTERWORTH;
};
ampConfig aFile;



// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
AmplitudeProcessor_Md::AmplitudeProcessor_Md() :
		AmplitudeProcessor("Md") {

	setSignalStart(0.);
	setSignalEnd(aFile.SIGNAL_LENGTH);
	setMinSNR(aFile.SNR_MIN);
	setMaxDist(8);
	_computeAbsMax = true;
	_isInitialized = false;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool AmplitudeProcessor_Md::setup(const Settings& settings) {

	if ( !AmplitudeProcessor::setup(settings) )
		return false;

	bool isButterworth = false;
	try {
		aFile.SEISMO = getParam<int>(settings, "md.seismo");
		std::string type;
		switch ( aFile.SEISMO ) {
			case 1:
				type = "WoodAnderson";
				break;
			case 2:
				type = "Seismo5sec";
				break;
			case 3:
				type = "WWSSN LP";
				break;
			case 4:
				type = "WWSSN SP";
				break;
			case 5:
				type = "Generic Seismometer";
				break;
			case 6:
				type = "Butterworth Low Pass";
				isButterworth = true;
				break;
			case 7:
				type = "Butterworth High Pass";
				isButterworth = true;
				break;
			case 8:
				type = "Butterworth Band Pass";
				isButterworth = true;
				break;
			case 9:
				type = "L4C 1Hz Seismometer";
				break;
			default:
				break;
		}
		SEISCOMP_DEBUG("md: %s sets SEISMO to %s [%s.%s]", AMPTAG, type.c_str(),
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.SEISMO = _SEISMO;
		SEISCOMP_ERROR("md: %s cannot read SEISMO value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	if ( isButterworth == true ) {
		try {
			aFile.BUTTERWORTH = getParam<std::string>(settings, "md.butterworth");
			SEISCOMP_DEBUG("md: %s sets Butterworth filter to  %s [%s.%s]", AMPTAG,
			               aFile.BUTTERWORTH.c_str(), settings.networkCode.c_str(),
			               settings.stationCode.c_str());
		}
		catch ( ... ) {
			aFile.BUTTERWORTH = _BUTTERWORTH;
			SEISCOMP_ERROR("md: %s cannot read Butterworth filter value from "
			               "configuration file [%s.%s]", AMPTAG,
			               settings.networkCode.c_str(), settings.stationCode.c_str());
		}
	}

	try {
		aFile.DEPTH_MAX = getParam<double>(settings, "md.depthmax");
		SEISCOMP_DEBUG("md: %s sets DEPTH MAX to  %.2f [%s.%s]", AMPTAG, aFile.DEPTH_MAX,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.DEPTH_MAX = _DEPTH_MAX;
		SEISCOMP_ERROR("md: %s cannot read DEPTH MAX value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		aFile.DELTA_MAX = getParam<double>(settings, "md.deltamax");
		SEISCOMP_DEBUG("md: %s sets DELTA MAX to  %.2f [%s.%s]", AMPTAG, aFile.DELTA_MAX,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.DELTA_MAX = _DELTA_MAX;
		SEISCOMP_ERROR("md: %s cannot read DELTA MAX value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		aFile.SNR_MIN = getParam<double>(settings, "md.snrmin");
		SEISCOMP_DEBUG("md: %s sets SNR MIN to  %.2f [%s.%s]", AMPTAG, aFile.SNR_MIN,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.SNR_MIN = _SNR_MIN;
		SEISCOMP_ERROR("md: %s cannot read SNR MIN value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		aFile.TAPER = getParam<double>(settings, "md.taper");
		SEISCOMP_DEBUG("md: %s sets TAPER to  %.2f [%s.%s]", AMPTAG, aFile.TAPER,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.TAPER = _TAPER;
		SEISCOMP_ERROR("md: %s cannot read TAPER value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		aFile.SIGNAL_LENGTH = getParam<double>(settings, "md.signal_length");
		SEISCOMP_DEBUG("md: %s sets SIGNAL LENGTH to  %.2f [%s.%s]", AMPTAG, aFile.SIGNAL_LENGTH,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.SIGNAL_LENGTH = _SIGNAL_LENGTH;
		SEISCOMP_ERROR("md: %s cannot read SIGNAL LENGTH value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		aFile.MD_MAX = getParam<double>(settings, "md.mdmax");
		SEISCOMP_DEBUG("md: %s sets MD MAX to  %.2f [%s.%s]", AMPTAG, aFile.MD_MAX,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.MD_MAX = _MD_MAX;
		SEISCOMP_ERROR("md: %s cannot read MD MAX value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		aFile.FMA = getParam<double>(settings, "md.fma");
		SEISCOMP_DEBUG("md: %s sets FMA to  %.4f [%s.%s]", AMPTAG, aFile.FMA,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.FMA = _FMA;
		SEISCOMP_ERROR("md: %s cannot read FMA value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		aFile.FMB = getParam<double>(settings, "md.fmb");
		SEISCOMP_DEBUG("md: %s sets FMB to  %.4f [%s.%s]", AMPTAG, aFile.FMB,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.FMB = _FMB;
		SEISCOMP_ERROR("md: %s cannot read FMB value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		aFile.FMD = getParam<double>(settings, "md.fmd");
		SEISCOMP_DEBUG("md: %s sets FMD to  %.4f [%s.%s]", AMPTAG, aFile.FMD,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.FMD = _FMD;
		SEISCOMP_ERROR("md: %s cannot read FMD value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		aFile.FMF = getParam<double>(settings, "md.fmf");
		SEISCOMP_DEBUG("md: %s sets FMF to  %.4f [%s.%s]", AMPTAG, aFile.FMF,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.FMF = _FMF;
		SEISCOMP_ERROR("md: %s cannot read FMF value from configuration file [%s.%s]", AMPTAG,
		               settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		aFile.FMZ = getParam<double>(settings, "md.fmz");
		SEISCOMP_DEBUG("md: %s sets FMZ to  %.4f [%s.%s]", AMPTAG, aFile.FMZ,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		aFile.FMZ = _FMZ;
		SEISCOMP_ERROR("md: %s cannot read FMZ value from configuration file [%s.%s]",
		               AMPTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	_isInitialized = true;

	return true;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void AmplitudeProcessor_Md::initFilter(double fsamp) {

	if ( !_enableResponses ) {

		SEISCOMP_DEBUG("Using custom responses");

		Math::Filtering::InPlaceFilter<double>* f;
		switch ( aFile.SEISMO ) {
			case 1:
				AmplitudeProcessor::setFilter(new Filtering::IIR::WoodAndersonFilter<
				        double>(Velocity, _config.woodAndersonResponse));
			break;
			case 2:
				AmplitudeProcessor::setFilter(new Filtering::IIR::Seismometer5secFilter<
				        double>(Velocity));
			break;
			case 3:
				AmplitudeProcessor::setFilter(new Filtering::IIR::WWSSN_LP_Filter<
				        double>(Velocity));
			break;
			case 4:
				AmplitudeProcessor::setFilter(new Filtering::IIR::WWSSN_SP_Filter<
				        double>(Velocity));
			break;
			case 5:
				AmplitudeProcessor::setFilter(new Filtering::IIR::GenericSeismometer<
				        double>(Velocity));
			break;
			case 6:
				f = new Math::Filtering::IIR::ButterworthLowpass<double>(3, 1, 15);
				AmplitudeProcessor::setFilter(f);
			break;
			case 7:
				f = new Math::Filtering::IIR::ButterworthHighpass<double>(3, 1, 15);
				AmplitudeProcessor::setFilter(f);
			break;
			case 8:
				// hardcoded ! We have to read the aFile.BUTTERWORTH
				f = new Math::Filtering::IIR::ButterworthBandpass<double>(3, 1, 15, 1);
				AmplitudeProcessor::setFilter(f);
			break;
			case 9:
				AmplitudeProcessor::setFilter(new Filtering::IIR::L4C_1Hz_Filter<
				        double>(Velocity));
			break;
			default:
				SEISCOMP_ERROR("md: %s cannot initialize the chosen filter, "
					"please review your configuration file", AMPTAG);
			break;
		}
	}
	else
		AmplitudeProcessor::setFilter(NULL);

	AmplitudeProcessor::initFilter(fsamp);
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int AmplitudeProcessor_Md::capabilities() const {
	return MeasureType;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
AmplitudeProcessor::IDList
AmplitudeProcessor_Md::capabilityParameters(Capability cap) const {

	if ( cap == MeasureType ) {
		IDList params;
		params.push_back("AbsMax");
		params.push_back("MinMax");
		return params;
	}

	return AmplitudeProcessor::capabilityParameters(cap);
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool AmplitudeProcessor_Md::setParameter(Capability cap,
                                         const std::string& value) {

	if ( cap == MeasureType ) {

		if ( value == "AbsMax" ) {
			_computeAbsMax = true;
			return true;
		}
		else if ( value == "MinMax" ) {
			_computeAbsMax = false;
			return true;
		}
		return false;
	}

	return AmplitudeProcessor::setParameter(cap, value);
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool AmplitudeProcessor_Md::deconvolveData(Response* resp,
                                           DoubleArray& data,
                                           int numberOfIntegrations) {
	if ( numberOfIntegrations < -1 )
		return false;

	SEISCOMP_DEBUG("Inside deconvolve function");

	double m, n;
	Math::Restitution::FFT::TransferFunctionPtr tf =
		resp->getTransferFunction(numberOfIntegrations < 0 ? 0 : numberOfIntegrations);

	if ( !tf )
		return false;

	Math::GroundMotion gm;

	if ( numberOfIntegrations < 0 )
		gm = Math::Displacement;
	else
		gm = Math::Velocity;

	Math::Restitution::FFT::TransferFunctionPtr cascade;
	Math::SeismometerResponse::WoodAnderson woodAndersonResp(gm, _config.woodAndersonResponse);
	Math::SeismometerResponse::Seismometer5sec seis5sResp(gm);
	Math::SeismometerResponse::L4C_1Hz l4c1hzResp(gm);

	Math::Restitution::FFT::PolesAndZeros woodAnderson(woodAndersonResp);
	Math::Restitution::FFT::PolesAndZeros seis5sec(seis5sResp);
	Math::Restitution::FFT::PolesAndZeros l4c1hz(l4c1hzResp);

	SEISCOMP_DEBUG("SEISMO = %d", aFile.SEISMO);

	switch ( aFile.SEISMO ) {
		case 1:
			cascade = *tf / woodAnderson;
		break;
		case 2:
			cascade = *tf / seis5sec;
		break;
		case 9:
			SEISCOMP_INFO("md: %s Applying filter L4C 1Hz to data", AMPTAG);
			cascade = *tf / l4c1hz;
		break;
		default:
			cascade = tf;
			SEISCOMP_INFO("md: %s No seismometer specified, no signal reconvolution performed", AMPTAG);
			return false;
		break;
	}

	// Remove linear trend
	Math::Statistics::computeLinearTrend(data.size(), data.typedData(), m, n);
	Math::Statistics::detrend(data.size(), data.typedData(), m, n);

    _config.respTaper = aFile.TAPER;
    SEISCOMP_DEBUG("md: %s TAPER is set to %.2f", AMPTAG, aFile.TAPER);

	return Math::Restitution::transformFFT(data.size(), data.typedData(),
	    _stream.fsamp, cascade.get(), _config.respTaper, _config.respMinFreq,
	    _config.respMaxFreq);
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool AmplitudeProcessor_Md::computeAmplitude(const DoubleArray& data, size_t i1,
                                             size_t i2, size_t si1, size_t si2,
                                             double offset, AmplitudeIndex* dt,
                                             AmplitudeValue* amplitude,
                                             double* period, double* snr) {

	double amax, Imax, ofs_sig, amp_sig = 0.0;
	DoubleArrayPtr d;

	if ( *snr < aFile.SNR_MIN )
		SEISCOMP_DEBUG("md: %s computed SNR is under configured SNR MIN", AMPTAG);

	if ( _computeAbsMax ) {
		size_t imax = find_absmax(data.size(), data.typedData(), si1, si2, offset);
		amax = fabs(data[imax] - offset);
		dt->index = imax;
	}
	else {
		int lmin, lmax;
		find_minmax(lmin, lmax, data.size(), data.typedData(), si1, si2, offset);
		amax = data[lmax] - data[lmin];
		dt->index = (lmin + lmax) * 0.5;
		dt->begin = lmin - dt->index;
		dt->end = lmax - dt->index;
	}

	Imax = dt->index;

	SEISCOMP_DEBUG("md: %s Amplitude max: %.2f", AMPTAG, amax);

	//! searching for Coda second by second through the end of the window
	//! if snrMin config is not 0 (config file or waveform review window)
	//! TODO: elevate accuracy by using a nanometers scale (maybe)
	if ( _config.snrMin != 0 ) {

		unsigned int i = si1;
		bool hasEndSignal = false;
		double calculatedSnr = -1;

		for ( i = (int) Imax; i < i2; i = i + 1 * (int) _stream.fsamp ) {

			int window_end = i + 1 * (int) _stream.fsamp;
			d = static_cast<DoubleArray*>(data.slice(i, window_end));

			//! computes pre-arrival offset
			ofs_sig = d->median();

			//! computes rms after removing offset
			amp_sig = 2 * d->rms(ofs_sig);

			if ( amp_sig / *_noiseAmplitude <= _config.snrMin ) {
				SEISCOMP_DEBUG("md: %s End of signal found! (%.2f <= %.2f)", AMPTAG,
				    (amp_sig / *_noiseAmplitude), _config.snrMin);
				hasEndSignal = true;
				calculatedSnr = amp_sig / *_noiseAmplitude;
				break;
			}
		}

		if ( !hasEndSignal ) {
			SEISCOMP_ERROR("md: %s SNR stayed over configured SNR_MIN! (%.2f > %.2f), "
				"skipping magnitude calculation for this station", AMPTAG,
			    calculatedSnr, _config.snrMin);
			return false;
		}

		dt->index = i;
	}
	else {
		dt->index = Imax;
	}

	//amplitude->value = 2 * amp_sig; //! actually it would have to be max. peak-to-peak
	amplitude->value = amp_sig;

	if ( _streamConfig[_usedComponent].gain != 0.0 ) {
		amplitude->value /= _streamConfig[_usedComponent].gain;
	}
	else {
		setStatus(MissingGain, 0.0);
		return false;
	}

	// Convert m/s to nm/s
	amplitude->value *= 1.E09;

	*period = dt->index - i1 + (_config.signalBegin * _stream.fsamp);

	SEISCOMP_DEBUG("md: %s calculated event amplitude = %.2f", AMPTAG, amplitude->value);
	SEISCOMP_DEBUG("md: %s calculated signal end at %.2f ms from P phase", AMPTAG, *period);

	return true;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
double AmplitudeProcessor_Md::timeWindowLength(double distance_deg) const {

	if ( !_isInitialized ) {

		aFile.MD_MAX = _MD_MAX;
		aFile.FMA = _FMA;
		aFile.FMZ = _FMZ;
		aFile.DEPTH_MAX = _DEPTH_MAX;
		aFile.STACOR = _STACOR;
		aFile.FMD = _FMD;
		aFile.FMB = _FMB;
		aFile.FMF = _FMF;
		aFile.SNR_MIN = _SNR_MIN;
		aFile.TAPER = _TAPER;
		aFile.SIGNAL_LENGTH = _SIGNAL_LENGTH;
		aFile.DELTA_MAX = _DELTA_MAX;
		aFile.SIGNAL_WINDOW_END = _SIGNAL_WINDOW_END;
		aFile.SEISMO = _SEISMO;
		aFile.BUTTERWORTH = _BUTTERWORTH;
	}

	double distance_km = Math::Geo::deg2km(distance_deg);
	double windowLength = (aFile.MD_MAX - aFile.FMA - (aFile.FMZ * aFile.DEPTH_MAX)
	        - aFile.STACOR - (aFile.FMD * distance_km)) / (aFile.FMB + aFile.FMF);

	windowLength = pow(10, windowLength) + aFile.SIGNAL_WINDOW_END;
	SEISCOMP_DEBUG("md: %s Requesting stream of %.2fsec for current station", AMPTAG, windowLength);

	return windowLength;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

/*----[ END OF AMPLITUDE PROCESSOR CLASS ]----*/




/*----[ MAGNITUDE PROCESSOR CLASS ]----*/

REGISTER_MAGNITUDEPROCESSOR(MagnitudeProcessor_Md, "Md");

struct magConfig {
	double DEPTH_MAX;
	double LINEAR_CORRECTION;
	double OFFSET;
	double DELTA_MAX;
	double MD_MAX;
	double TAPER;
	double SIGNAL_LENGTH;
	double SNR_MIN;
	double FMA;
	double FMB;
	double FMD;
	double FMF;
	double FMZ;
	double STACOR;
};

magConfig mFile;



// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
MagnitudeProcessor_Md::MagnitudeProcessor_Md() :
		MagnitudeProcessor("Md") {

	_linearCorrection = mFile.LINEAR_CORRECTION;
	_constantCorrection = mFile.OFFSET;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool MagnitudeProcessor_Md::setup(const Settings& settings) {

	try {
		mFile.DELTA_MAX = getParam<double>(settings, "md.deltamax");
		SEISCOMP_DEBUG("md: %s sets DELTA MAX to  %.2f [%s.%s]", MAGTAG, mFile.DELTA_MAX,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.DELTA_MAX = _DELTA_MAX;
		SEISCOMP_ERROR("md: %s cannot read DELTA MAX value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.DEPTH_MAX = getParam<double>(settings, "md.depthmax");
		SEISCOMP_DEBUG("md: %s sets DEPTH MAX to  %.2f [%s.%s]", MAGTAG, mFile.DEPTH_MAX,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.DEPTH_MAX = _DEPTH_MAX;
		SEISCOMP_ERROR("md: %s cannot read DEPTH MAX value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.MD_MAX = getParam<double>(settings, "md.mdmax");
		SEISCOMP_DEBUG("md: %s sets MD MAX to  %.2f [%s.%s]", MAGTAG, mFile.MD_MAX,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.MD_MAX = _MD_MAX;
		SEISCOMP_ERROR("md: %s cannot read MD MAX value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.LINEAR_CORRECTION = getParam<double>(settings, "md.linearcorrection");
		SEISCOMP_DEBUG("md: %s sets LINEAR CORRECTION to  %.2f [%s.%s]", MAGTAG,
		    mFile.LINEAR_CORRECTION, settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.LINEAR_CORRECTION = _LINEAR_CORRECTION;
		SEISCOMP_ERROR("md: %s cannot read LINEAR CORRECTION value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.OFFSET = getParam<double>(settings, "md.offset");
		SEISCOMP_DEBUG("md: %s sets OFFSET to  %.2f [%s.%s]", MAGTAG, mFile.OFFSET,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.OFFSET = _OFFSET;
		SEISCOMP_ERROR("md: %s cannot read OFFSET value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.FMA = getParam<double>(settings, "md.fma");
		SEISCOMP_DEBUG("md: %s sets FMA to  %.4f [%s.%s]", MAGTAG, mFile.FMA,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.FMA = _FMA;
		SEISCOMP_ERROR("md: %s cannot read FMA value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.FMB = settings.getDouble("magnitudes.md.fmb");
		SEISCOMP_DEBUG("md: %s sets FMB to  %.4f [%s.%s]", MAGTAG, mFile.FMB,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.FMB = _FMB;
		SEISCOMP_ERROR("md: %s cannot read FMB value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.FMD = settings.getDouble("magnitudes.md.fmd");
		SEISCOMP_DEBUG("md: %s sets FMD to  %.4f [%s.%s]", MAGTAG, mFile.FMD,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.FMD = _FMD;
		SEISCOMP_ERROR("md: %s cannot read FMD value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.FMF = settings.getDouble("magnitudes.md.fmf");
		SEISCOMP_DEBUG("md: %s sets FMF to  %.4f [%s.%s]", MAGTAG, mFile.FMF,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.FMF = _FMF;
		SEISCOMP_ERROR("md: %s cannot read FMF value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.SNR_MIN = getParam<double>(settings, "md.snrmin");
		SEISCOMP_DEBUG("md: %s sets SNR MIN to  %.4f [%s.%s]", MAGTAG, mFile.SNR_MIN,
				settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.SNR_MIN = _SNR_MIN;
		SEISCOMP_ERROR("md: %s cannot read SNR MIN value from configuration file [%s.%s]",
				MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	try {
		mFile.TAPER = getParam<double>(settings, "md.taper");
		SEISCOMP_DEBUG("md: %s sets TAPER to  %.4f [%s.%s]", MAGTAG, mFile.TAPER,
				settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.TAPER = _TAPER;
		SEISCOMP_ERROR("md: %s cannot read TAPER value from configuration file [%s.%s]",
				MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	try {
		mFile.SIGNAL_LENGTH = getParam<double>(settings, "md.signal_length");
		SEISCOMP_DEBUG("md: %s sets SIGNAL LENGTH to  %.4f [%s.%s]", MAGTAG, mFile.SIGNAL_LENGTH,
				settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.SIGNAL_LENGTH = _SIGNAL_LENGTH;
		SEISCOMP_ERROR("md: %s cannot read SIGNAL LENGTH value from configuration file [%s.%s]",
				MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.FMZ = getParam<double>(settings, "md.fmz");
		SEISCOMP_DEBUG("md: %s sets FMZ to  %.4f [%s.%s]", MAGTAG, mFile.FMZ,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.FMZ = _FMZ;
		SEISCOMP_ERROR("md: %s cannot read FMZ value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	try {
		mFile.STACOR = getParam<double>(settings, "md.stacor");
		SEISCOMP_DEBUG("md: %s sets STACOR to  %.4f [%s.%s]", MAGTAG, mFile.STACOR,
		    settings.networkCode.c_str(), settings.stationCode.c_str());
	}
	catch ( ... ) {
		mFile.STACOR = _STACOR;
		SEISCOMP_ERROR("md: %s cannot read STACOR value from configuration file [%s.%s]",
		    MAGTAG, settings.networkCode.c_str(), settings.stationCode.c_str());
	}

	return true;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
MagnitudeProcessor::Status
MagnitudeProcessor_Md::computeMagnitude(double amplitude, const std::string &,
                                        double period, double, double delta, double depth,
                                        const DataModel::Origin *,
                                        const DataModel::SensorLocation *,
                                        const DataModel::Amplitude *,
                                        const Locale *,
                                        double& value) {

	double epdistkm;
	epdistkm = Math::Geo::deg2km(delta);

	SEISCOMP_DEBUG("md: %s --------------------------------", MAGTAG);
	SEISCOMP_DEBUG("md: %s |    PARAMETERS   |    VALUE   |", MAGTAG);
	SEISCOMP_DEBUG("md: %s --------------------------------", MAGTAG);
	SEISCOMP_DEBUG("md: %s | window length   | %.2f ", MAGTAG, mFile.SIGNAL_LENGTH);
	SEISCOMP_DEBUG("md: %s | taper           | %.2f ", MAGTAG, mFile.TAPER);
	SEISCOMP_DEBUG("md: %s | min snr         | %.2f ", MAGTAG, mFile.SNR_MIN);
	SEISCOMP_DEBUG("md: %s | delta max       | %.2f ", MAGTAG, mFile.DELTA_MAX);
	SEISCOMP_DEBUG("md: %s | depth max       | %.2f ", MAGTAG, mFile.DEPTH_MAX);
	SEISCOMP_DEBUG("md: %s | md max          | %.2f ", MAGTAG, mFile.MD_MAX);
	SEISCOMP_DEBUG("md: %s | fma             | %.4f ", MAGTAG, mFile.FMA);
	SEISCOMP_DEBUG("md: %s | fmb             | %.4f ", MAGTAG, mFile.FMB);
	SEISCOMP_DEBUG("md: %s | fmd             | %.4f ", MAGTAG, mFile.FMD);
	SEISCOMP_DEBUG("md: %s | fmf             | %.4f ", MAGTAG, mFile.FMF);
	SEISCOMP_DEBUG("md: %s | fmz             | %.4f ", MAGTAG, mFile.FMZ);
	SEISCOMP_DEBUG("md: %s | stacor          | %.4f ", MAGTAG, mFile.STACOR);
	SEISCOMP_DEBUG("md: %s --------------------------------", MAGTAG);
	SEISCOMP_DEBUG("md: %s | (f-p)           | %.2f sec ", MAGTAG, period);
	SEISCOMP_DEBUG("md: %s | seismic depth   | %.2f km ", MAGTAG, depth);
	SEISCOMP_DEBUG("md: %s | epicenter dist  | %.2f km ", MAGTAG, epdistkm);
	SEISCOMP_DEBUG("md: %s --------------------------------", MAGTAG);

	if ( amplitude <= 0. ) {
		value = 0;
		SEISCOMP_ERROR("md: %s calculated amplitude is wrong, "
		               "no magnitude will be calculated", MAGTAG);
		return Error;
	}

	if ( (mFile.DELTA_MAX) < epdistkm ) {
		SEISCOMP_ERROR("md: %s epicenter distance is out of configured range, "
		               "no magnitude will be calculated", MAGTAG);
		return DistanceOutOfRange;
	}

	if ( period <= 0 ) {
		return PeriodOutOfRange;
	}

	value = mFile.FMA + mFile.FMB * log10(period) + (mFile.FMF * period)
	        + (mFile.FMD * epdistkm) + (mFile.FMZ * depth) + mFile.STACOR;

	if ( value > mFile.MD_MAX )
		SEISCOMP_WARNING("md: %s Calculated magnitude is beyond max Md value [value= %.2f]",
		                 MAGTAG, value);

	return OK;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

/*----[ END OF MAGNITUDE PROCESSOR CLASS ]----*/


