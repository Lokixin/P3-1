/// @file

#include <iostream>
#include <fstream>
#include <string.h>
#include <errno.h>
#include <bits/stdc++.h> 
#include "wavfile_mono.h"
#include "pitch_analyzer.h"

#include "docopt.h"

#define FRAME_LEN   0.030 /* 30 ms. */
#define FRAME_SHIFT 0.015 /* 15 ms. */

using namespace std;
using namespace upc;

static const char USAGE[] = R"(
get_pitch - Pitch Detector 

Usage:
    get_pitch [options] <input-wav> <output-txt>
    get_pitch (-h | --help)
    get_pitch --version

Options:
    -h, --help  Show this screen
    --version   Show the version of the project

Arguments:
    input-wav   Wave file with the audio signal
    output-txt  Output file: ASCII file with the result of the detection:
                    - One line per frame with the estimated f0
                    - If considered unvoiced, f0 must be set to f0 = 0
)";

int main(int argc, const char *argv[]) {
	/// \TODO 
	///  Modify the program syntax and the call to **docopt()** in order to
	///  add options and arguments to the program.
    std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
        {argv + 1, argv + argc},	// array of arguments, without the program name
        true,    // show help if requested
        "2.0");  // version string

	std::string input_wav = args["<input-wav>"].asString();
	std::string output_txt = args["<output-txt>"].asString();

  // Read input sound file
  unsigned int rate;
  vector<float> x;
  if (readwav_mono(input_wav, rate, x) != 0) {
    cerr << "Error reading input file " << input_wav << " (" << strerror(errno) << ")\n";
    return -2;
  }

  int n_len = rate * FRAME_LEN;
  int n_shift = rate * FRAME_SHIFT;

  // Define analyzer
  PitchAnalyzer analyzer(n_len, rate, PitchAnalyzer::HAMMING, 50, 500);

  /// \TODO
  /// Preprocess the input signal in order to ease pitch estimation. For instance,
  /// central-clipping or low pass filtering may be used.

    //Low Pass Filter
  /*float alph1 = 0.7; 
  for(unsigned int low= 2; low < x.size(); ++low){
    x[low]= alph1*x[low] + (1-alph1)*x[low-1];
  }*/
  int M = 7;
  float aux = 0;
  for(int i = M/2; i<x.size() - M/2; i++){
    for(int j = -M/2; j<=M/2; j++){
      aux += x[i+j]/M;
      x[i] = aux;
    }aux = 0;
  }

  //Central clipping
  unsigned int WINDOW= 200;
  double maximo = *max_element(x.begin(), x.end());
  for(unsigned int clip = 0; clip < x.size(); ++clip){
    if(0.05*maximo > x[clip]) x[clip] = 0;
  }
  for(unsigned int idx = 0; idx < x.size()/WINDOW - 1; ++idx){
    for (unsigned int w_idx = 0; idx < WINDOW; ++idx){
      maximo = *max_element(x.begin() + idx*WINDOW, x.begin() + (idx+1)*WINDOW);
      if(maximo*0.15>x[w_idx + idx*WINDOW]){
        x[w_idx + idx*WINDOW] = 0;
      }
      
    }
  }

  // Iterate for each frame and save values in f0 vector
  vector<float>::iterator iX;
  vector<float> f0;
  for (iX = x.begin(); iX + n_len < x.end(); iX = iX + n_shift) {
    float f = analyzer(iX, iX + n_len);
    f0.push_back(f);
  }


  /// \TODO
  /// Postprocess the estimation in order to supress errors. For instance, a median filter
  /// or time-warping may be used.

  //median filter
  int median_size = 5;
  int median_center = median_size/2;
  float m_w[median_size];
  for(unsigned int pos = median_center; pos < f0.size()-median_center; ++pos){
    for(int med_idx = -median_center; med_idx <= median_center; ++med_idx){
      m_w[med_idx + median_center] = f0[pos + med_idx];
    }
    sort(m_w, m_w+median_size);
    f0[pos] = m_w[median_center];
  }

  // Write f0 contour into the output file
  ofstream os(output_txt);
  if (!os.good()) {
    cerr << "Error reading output file " << output_txt << " (" << strerror(errno) << ")\n";
    return -3;
  }

  os << 0 << '\n'; //pitch at t=0
  for (iX = f0.begin(); iX != f0.end(); ++iX) 
    os << *iX << '\n';
  os << 0 << '\n';//pitch at t=Dur

  return 0;
}
