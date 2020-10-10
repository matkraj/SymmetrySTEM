#ifndef SYMMETRY_SUPPORT_H
#define SYMMETRY_SUPPORT_H

#include <cstdio>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <time.h>       /* time_t, struct tm, time, localtime, strftime */
#include <fstream>
#include <algorithm>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <arrayfire.h>

#include <dirent.h>
#include <cstring>
#include <random>


struct DefaultVals{
    std::string dirname,binfile;
    bool DOG,endian,sign,show,AUTO,hann,gradient;	 	
    unsigned int HEADERSIZE,FOOTERSIZE;	     
    int STEP,ind;		
    unsigned long ZEROPAD, SIZE_DET, SCAN_X, SCAN_Y, SCAN_XY, BYTESIZE_DET, frame_i,leftovers;	
    double sigma,dogSmall,dogLarge,DOWNSCALE,angle;
    uint bitdepth,bytedepth;
    af::timer start;
};

struct Results{
    af::array sums,rotsim,rotmirror,pacbed;
};

void showhelpinfo();
void showMe();
void gotHere();
void calcString();
float roundFloat(float in);
void data_check(std::ifstream& fs, DefaultVals s);
af::array normalizeImage(af::array in);
af::array normalizeFloat(af::array in);
void showImage(af::array in, af::Window& myWindow);
af::array pyramid(const af::array& img, const int level);
bool DirectoryExists(const char* pzPath );
std::string sConvert (float number);
af::array readDiskLine(std::ifstream& fs, DefaultVals s);


int compute(std::ifstream& fs, DefaultVals &s, Results &calcs, af::Window& myWindow, bool end);

void generate_outputs(DefaultVals s,Results &calcs);
Results createResults(DefaultVals s);
void createLogFile(DefaultVals s);
int getSettings(DefaultVals &s, int argc, char** argv);
af::array generate_edge(DefaultVals s, af::array mask);
af::array generate_edge(af::array disk);
af::array partialConvolve(DefaultVals s, af::array edge_line,af::array EdgeLong);
af::array cross_correlate(const af::array &a1, const af::array &a2);
af::array normalizeCC(const af::array &in);
af::array smooth(const af::array &in, float smooth = 1.0, af::array smooths = af::constant(1,1));
af::array smooth1D(const af::array &in, float smooth = 1.0);
void normalizeCCinPlace(af::array &in);
void normaliseViewInPlace(af::array& disk);
af::array normaliseView(const af::array& disk);
af::array generate_edge2(const af::array &disk, float smooth = 1.0);
af::array zoom_image(const af::array &disk, float zoom);
af::array maxnorm(const af::array &disk);
af::array get_disks_RANDOM(std::ifstream& fs,DefaultVals s, unsigned long noFrames);
af::array cross_correlate2(const af::array &a1, const af::array &a2);
af::array generate_hann(const af::array &disk);
void plot1d(af::Window &win, af::array ar1d, std::string lala = "1d data");
af::array removeHotPixels(const af::array &in, float thr);
af::array zeropad(const DefaultVals &s, const af::array &ar_line_small);
af::array zeropad(float scale, const af::array &ar_line_small);
af::array normalize3D(const af::array &in);
void endianSwap(DefaultVals s, int LENGTH, char* tmp);
af::array pointerToArray(DefaultVals s, char* tmpNoHead, bool mask = false);

    
void check_s(DefaultVals &s);
af::array get_mask(DefaultVals s);
void setDefaultSettings(DefaultVals &s);
void debug_win(af::array a,af::array b, af::array c, af::array d,af::array e,af::array f,af::array g,af::array h, af::Window& myWindow);
void saveimages(DefaultVals s, Results calcs);
af::array norm_hist(af::array in, float percentage);
af::array smoothed_disk(int size, int radius, float smooth);

#endif
