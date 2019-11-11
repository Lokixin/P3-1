/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"
#include <ffft/FFTReal.h>
#include <iomanip>
using namespace std;

/// Name space of UPC
namespace upc {
/// \TODO Compute the autocorrelation r[l]
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {
    for (unsigned int k = 0; k < r.size(); ++k) {
      r[k] = 0;
      for (unsigned int i = 0; i< r.size()-k; ++i){
        r[k] += x[i]*x[i+k];
      }
      r[k] = r[k]/(r.size()-k);
    }

    if (r[0] == 0.0F) //to avoid log() and divide zero 
      r[0] = 1e-10; 
  }
  void PitchAnalyzer::amdf(const vector<float> &x, vector<float> &s) const {
      for (unsigned int k = 0; k < s.size(); ++k) {
        for (unsigned int i = 0; i< s.size()-k; ++i){
          s[k] += abs(x[i]-x[i+k]);
        }
      }
    }


  void PitchAnalyzer::cepstrum(const vector<float> &x, vector<float> &c) const {
    
    unsigned int N = 0;
    unsigned int i = 2;
    
    
    while (N < x.size()){
      N = pow(N,i);
      i += 1;
    }
    ffft::FFTReal <float> fft_object(N);
    float aux[N]; 
    for(unsigned int n =0; n < N; ++n){
      if(n < x.size()){
        aux[n] = x[n];
      }else{
        aux[n] = 0;
      }
    }
    float out[N];
    fft_object.do_fft (out, aux);  
    for (unsigned int m = 0; m < N; ++m){
      out[m] = log10(abs(out[m]));
    }
    float in[N];
    fft_object.do_ifft (out, in);   
    fft_object.rescale (in); 
    for (unsigned int o = 0; o < N; ++o){
      c[o] = in[o];
    }
  }

  void PitchAnalyzer::set_window(Window win_type) {
    if (frameLen == 0)
      return;

    window.resize(frameLen);

    switch (win_type) {
    case HAMMING:
      /// \TODO Implement the Hamming window
      for (unsigned int i = 0; i < window.size(); ++i){
        window[i] = 0.54 - 0.46*cos(2*M_PI*i/frameLen);
      }
      break;
    case RECT:
    default:
      window.assign(frameLen, 1);
    }
  }

  void PitchAnalyzer::set_f0_range(float min_F0, float max_F0) {
    npitch_min = (unsigned int) samplingFreq/max_F0;
    if (npitch_min < 2)
      npitch_min = 2;  // samplingFreq/2

    npitch_max = 1 + (unsigned int) samplingFreq/min_F0;

    //frameLen should include at least 2*T0
    if (npitch_max > frameLen/2)
      npitch_max = frameLen/2;
  }
  float PitchAnalyzer::zeros(vector<float> & x)const {
    float zeros = 0;
    for(unsigned int i = 1; i < x.size(); ++i){
      if((x[i-1] <= 0 && x[i] > 0) || (x[i-1] >= 0 && x[i] < 0)){
        zeros +=1;
      }    
    }
    
    return zeros;
  }

  bool PitchAnalyzer::unvoiced(float pot, float r1norm, float rmaxnorm, float zeros) const {
    /// \TODO Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones. 
    //false significa sonora, true significa sorda
    float th_1 = 0.95;
    float th_2 = 0.63;
    float th_zcr = 1100;
    float th_pot = -48;
    if((r1norm >= th_1 &&  pot >= th_pot) || (rmaxnorm >= th_2 && zeros <= th_zcr)){
      return false;
    }
    else {
      return true;
    }
    /*if(r1norm >= 0.82  && pot >= -65) return false;
    else return true;*/
    
  }

  float PitchAnalyzer::compute_pitch(vector<float> & x) const {
    if (x.size() != frameLen)
      return -1.0F;

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i)
      x[i] *= window[i];

    vector<float> r(npitch_max);
    vector<float> s(npitch_max);
    //Compute correlation
    autocorrelation(x, r);
    amdf(x,s);
    //cepstrum(x, c);
    
    vector<float>::const_iterator iR = r.begin(), iRMax = iR;

    /// \TODO 
	/// Find the lag of the maximum value of the autocorrelation away from the origin.<br>
	/// Choices to set the minimum value of the lag are:
	///    - The first negative value of the autocorrelation.
	///    - The lag corresponding to the maximum value of the pitch.
    ///	   .
	/// In either case, the lag should not exceed that of the minimum value of the pitch.

    unsigned int lag = iRMax - r.begin();
    unsigned int pmin = 60;
    unsigned int pmax = 300;
    
    lag = max_element(r.begin()+pmin, r.begin()+pmax) - r.begin();
    //lag = min_element(s.begin()+pmin, s.begin()+pmax) - s.begin();
    float pot = 10 * log10(r[0]);

    //You can print these (and other) features, look at them using wavesurfer
    //Based on that, implement a rule for unvoiced
    //change to #if 1 and compile
    double zero = samplingFreq*zeros(x)/(2*(frameLen-1));

#if 1
    if (r[0] > 0.0F){
      cout << setprecision(2)<< zero <<'\t'<< r[1]/r[0] <<'\t'<< r[lag]/r[0] << endl;
    }else{
      cout <<setprecision(2) << 0 <<'\t'<< 0 <<'\t'<< 0 << endl;
    }
      
#endif
    
    if (unvoiced(pot, r[1]/r[0], r[lag]/r[0], zero) )
      return 0;
    else
      return (float) samplingFreq/(float)lag;
  }
}
