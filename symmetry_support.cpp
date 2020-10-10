#include "symmetry_support.h"

using namespace std;
using namespace af;


// help output with a couple of examples
void showhelpinfo(){
  cout<<"Usage:   "<<" [-option] [argument]"<<endl;
  cout<<"option:  "<<"-h  show help information"<<endl;
  cout<<"         "<<"-o  output dir"<<endl;
  cout<<"         "<<"-b  binary big file 8/16/32 bit"<<endl;
  cout<<"         "<<"-A  angle to compute symmetry analysis"<<endl;
  cout<<"         "<<"=======DETECTOR/FILE SETTINGS======"<<endl;
  cout<<"         "<<"-e  size of header of detector frame in bytes in front of each frame (default 384)"<<endl;
  cout<<"         "<<"-c  size of footer of detector frame in bytes in back of each frame (default 0)"<<endl;  
  cout<<"         "<<"-p  number of frames processed in parallel (default 128)"<<endl;
  cout<<"         "<<"-x  'x' dimension of the image (default 257)"<<endl;
  cout<<"         "<<"-y  'y' dimension of the image (default 256)"<<endl;
  cout<<"         "<<"-z  zero padding (multiples of detector size)"<<endl;
  cout<<"         "<<"-d  size of square detector (default 256)"<<endl;
  cout<<"         "<<"-l  unsigned/signed switch - unsigned 0, signed 1 (default 1)"<<endl;
  cout<<endl<<"Example usage\t $ ./symmetry -b bin_file -o output_dir -e 0 -c 1024 -A 180 -d 128 -x 256 -y 256 -a 32 "
              "\n\n This computes a rotation and mirror symmetry analysis for 180 degree rotation angle for a RAW binary data of square "
              "detector with 128 pixels with 32 bitdepth, no header before each frame and 1024 bytes footer size.\n\n"<<endl;
}

float roundFloat(float in) {
    long inn = long(in*100);
    return inn/100.0;
}

af::array normalizeImage(array in) {
    float min = af::min<float>(in);
    float max = af::max<float>(in);
    return 65635.0f*((in - min) / (max - min));
}

af::array normalizeCC(const af::array &in) {
    unsigned int x = in.dims(0);
    unsigned int y = in.dims(1);
    unsigned int z = in.dims(2);
    // make2D array
    array in2d = af::moddims(in,x*y,z);
    array mean2d = af::mean(in2d,0);
    array stdev2d = af::stdev(in2d,0);
    in2d = ( in2d - af::tile(mean2d,x*y) ) / af::tile(stdev2d,x*y);
    return af::moddims(in2d,x,y,z);
}

af::array normalize3D(const af::array &in) {
    unsigned int x = in.dims(0);
    unsigned int y = in.dims(1);
    unsigned int z = in.dims(2);
    // make2D array
    
    array in2d = af::moddims(in,x*y,z);
    array min2d = af::min(in2d,0);
    in2d = in2d - af::tile(min2d,x*y,1);
    array max2d = af::max(in2d,0);
    in2d = in2d / af::tile(max2d,x*y,1);
    return af::moddims(in2d,x,y,z);
}

void normalizeCCinPlace(af::array &in) {
    unsigned int x = in.dims(0);
    unsigned int y = in.dims(1);
    unsigned int z = in.dims(2);

    in = af::moddims(in,x*y,z);
    array mean2d = af::mean(in,0);
    array stdev2d = af::stdev(in,0);
    

    in = ( in - af::tile(mean2d,x*y, 1) ) / af::tile(stdev2d,x*y,1);
    in = af::moddims(in,x,y,z);
}

af::array normalizeFloat(array in) {
    in = in.as(f32);
    float min = af::min<float>(in);
    float max = af::max<float>(in);
    in = (in - min) / (max - min);
    return in;
}

af::array generate_hann(const af::array &disk) {
    //this whole thing is probably wrong for 2D images
    af::array rangeX = af::range(disk.dims(0));
    af::array rangeY = af::range(disk.dims(1));

    rangeX = 0.5*(1-af::cos(2*M_PI*maxnorm(rangeX)));
    rangeY = 0.5*(1-af::cos(2*M_PI*maxnorm(rangeY)));
    
    rangeX = af::tile(rangeX,1,disk.dims(0));
    rangeY = af::reorder(rangeY,1,0);
    rangeY = af::tile(rangeY,disk.dims(1),1);
    
    af::array hann = rangeX * rangeY;
    hann = 1.0;
    return maxnorm(hann);
}

af::array norm_hist(array in, float percentage) {
    float mn = min<float>(in);
    float mx = max<float>(in);
    
    int noBins = 1024;
    array hist = histogram(in,noBins,mn,mx);
    
    int index = 0;
    int sumd = 0;
    float maxsum = in.dims(0) * in.dims(1) / 100 * percentage; // break for histogram stretch - to cover for errors in analysis
    
    //calculate minimum for resulting array
    while (sumd < maxsum) {
        sumd += sum<int>(hist(seq(index,index)));
        index++;
    }
    float minb = mn + index * ( (mx-mn) / noBins);

    index = 0;
    sumd = 0;
    //calculate maximum for resulting array
    while (sumd < maxsum) {
        sumd += sum<int>(hist(seq(noBins-index-1,noBins-index-1)));
        index++;
    }
    float maxb = mx - index * ( (mx-mn) / noBins);
    in = 65635.0f * ((in - minb) / (maxb - minb));
    return in;
}

void showImage(af::array in, af::Window &myWin) {
    while (not myWin.close()) {
        myWin.image(normaliseView(in));
    }
}

//function to upsample - may be used for other types of coding subpixel registration
array pyramid(const array& img, const int level)
{
    array pyr = img.copy();
    for(int i = 0; i < level; i++) {
        array tmp = constant(0, pyr.dims(0) * 2, pyr.dims(1) * 2, pyr.dims(2));
        tmp(seq(0, 2*pyr.dims(0)-1, 2), seq(0, 2*pyr.dims(1)-1, 2), span) = pyr;
        pyr=tmp;
    }
    return pyr;
}

af::array smoothed_disk(int size, int radius, float smooth) {
    if (radius == 0) {
        return af::gaussianKernel( size, size , smooth, smooth);
    } else {
        af::array x = range(size) - size/2;
        af::array y = af::reorder(x,1,0);
        
        x = af::tile(x,1,size);
        y = af::tile(y,size,1);
            
        af::array disk = (sqrt(x*x+y*y)<radius)*1.0f;
        if (smooth > 0.001) {
            af::array gaussian = af::gaussianKernel( size, size , smooth, smooth);
            disk = fftConvolve2(gaussian,disk,AF_CONV_DEFAULT);
        }
        return maxnorm(disk);
    }
}

void normaliseViewInPlace(af::array& disk) {
    disk = disk - af::min<float>(disk);
    disk = disk / af::max<float>(disk);
}

af::array normaliseView(const af::array& disk) {
    array disk2 = disk - af::min<float>(disk);
    disk2 = disk2 / af::max<float>(disk2);
    return disk2;
}

af::array normaliseViewRange1D(const af::array& in, unsigned long maxrange) {
    array out = in - af::min<float>(in(seq(0,maxrange)));
    out = out / af::max<float>(out(seq(0,maxrange)));
    return out;
}

af::array smooth(const af::array &in, float smooth, af::array smooths) {
    af::array gaussian;
    if (smooths.dims(0) == 1) {
        gaussian = af::gaussianKernel(in.dims(0),in.dims(1), smooth, smooth);
        return af::fftConvolve2(in,gaussian,AF_CONV_DEFAULT);
    } else {
        array smooths3D = af::constant(0,in.dims(0),in.dims(1),smooths.dims(0));
        for (int x=0;x<smooths.dims(0);x++) {
            gaussian = af::gaussianKernel(
                in.dims(0),in.dims(1), 
                af::sum<float>(smooths(seq(x,x))), af::sum<float>(smooths(seq(x,x)))
            );
            smooths3D.slice(x) = af::fftConvolve2(in,gaussian,AF_CONV_DEFAULT);
        }
        return smooths3D;
    }
}
    
af::array smooth1D(const af::array &in, float smooth) {
    af::array gaussian = af::gaussianKernel(in.dims(0),in.dims(0), smooth, smooth);
    gaussian = gaussian(seq(),seq(gaussian.dims(0)/2,gaussian.dims(0)/2));
    return af::fftConvolve2(in,gaussian,AF_CONV_DEFAULT );
}

bool DirectoryExists(const char* pzPath) {
    if ( pzPath == NULL) return false;
    DIR *pDir;
    bool bExists = false;
    pDir = opendir (pzPath);
    if (pDir != NULL){
        bExists = true;
        (void) closedir (pDir);
    }
    return bExists;
}

af::array maxnorm(const af::array &disk) {
    return disk.as(f32)/af::max<float>(disk);
}

unsigned long posMax1D(const array& in) {
    af::array maximum, position;
    af::max(maximum,position,af::flat(in));
    return af::sum<unsigned long>(position);
}

unsigned long posMin1D(const array& in) {
    af::array maximum, position;
    af::max(maximum,position,af::flat(in));
    return af::sum<unsigned long>(position);
}

af::array integral1D(const array& in) {
    af::array out = in.copy();
    for (int x=0;x<in.dims(0);x++)
        out(seq(x,x)) = af::sum<float>(in(seq(0,x)));
    return out;
}

af::array removeHotPixels(const af::array &in, float thr) {
    //compute histogram
    af::array hist = af::histogram(in,256);
    // remove hotest pixels 0.01% -ish
    // CUDA code faills if the following or a hanning window is not used
    float sumhist = af::sum<float>(hist);
    int r = hist.dims(0)-1;
    af::array subhist;
    //get rid of 0.1 % of hot pixels
    while (2>1) {
        subhist = hist(seq(r--,hist.dims(0)-1));
        float sumsub = af::sum<float>(subhist);
        if (sumsub/sumhist > thr) break;
    }
    
    float min_in = af::min<float>(in);
    float max_in = af::max<float>(in);
    float threshold = (max_in - min_in)*r/hist.dims(0) + min_in;
    
    af::array out = in.copy();
    out = (in<threshold)*in + (in>=threshold)*threshold;
    af::deviceGC();
    return out;
}

void plot1d(af::Window &win, af::array ar1d, std::string lala) {
    //plot 1D line graph
    af::array x = af::range(ar1d.dims(0));
    win.plot(x.as(f32),ar1d.as(f32),lala.c_str()); //scatter with AF_MARKER_CIRCLE
}

af::array cross_correlate2(const af::array &a1, const af::array &a2) { 
    //cross correlation with upscale
    unsigned int scale = 2;
    array A1 = af::fft2Norm(a1,double(scale),a1.dims(0)*scale,a1.dims(1)*scale);
    array A2 = af::fft2Norm(a2,double(scale),a2.dims(0)*scale,a2.dims(1)*scale);
    array ccor_line = af::real(af::ifft2(A1*af::conjg(A2)));
	// INVERT QUADRANTS AFTER IFFT
	ccor_line = shift(ccor_line,ccor_line.dims(0)/2,ccor_line.dims(1)/2);
    ccor_line /= float(ccor_line.dims(0))*float(ccor_line.dims(1));
    //ccor_line = ccor_line * scale * scale;
    unsigned int l = a1.dims(0)/2*(scale - 1);
    unsigned int h = a1.dims(0)/2*(scale + 1) -1;
    return ccor_line(seq(l,h),seq(l,h));
}


af::array get_disks_RANDOM(std::ifstream& fs, DefaultVals s,  unsigned long noFrames) {
    //calculate size of read out data for the *tmp0 pointer in bytes

    int FRAMELENGTH = s.HEADERSIZE+s.BYTESIZE_DET+s.FOOTERSIZE;
    int LENGTH = FRAMELENGTH * noFrames; 
    int DATALENGTH = s.BYTESIZE_DET * noFrames;

    //stack overflow answer to random number generation
    std::mt19937 rng;
    rng.seed(std::random_device()());
    // distribution to search for random set of frames not over the edge of acquisition
    
    unsigned long limit;
    unsigned long newpos;
    
    //this will bias random search in 'x' direction - however due to microscope stability and 
    //scanning it may not be a bad thing 
    // (also would not want to code reader in 'y' direction - hopping in wrong way )

    
    if (noFrames>=long(0.75*s.SCAN_X)) {
        limit = s.SCAN_XY;
        std::uniform_int_distribution<std::mt19937::result_type> dist6(0,limit - noFrames); 
    
        unsigned long r = dist6(rng);
        //somehow the random number overflows 
        while (r>s.SCAN_XY-noFrames) r = dist6(rng);

        newpos = FRAMELENGTH * r;
        
    } else { 
        limit = s.SCAN_X;
        if (limit>s.SCAN_Y) limit = s.SCAN_Y;
        std::uniform_int_distribution<std::mt19937::result_type> dist6(0,limit - noFrames); 
    
        unsigned long rx = dist6(rng);
        unsigned long ry = dist6(rng);
        //somehow the random number overflows 
        while (rx>s.SCAN_XY-noFrames) rx = dist6(rng);
        while (ry>s.SCAN_XY-noFrames) ry = dist6(rng);

        newpos = FRAMELENGTH * rx + FRAMELENGTH * ry * s.SCAN_X;
    }
    
    fs.seekg(newpos, ios::beg);

    char *tmp0 = new char [LENGTH];
    fs.read((char*)tmp0, LENGTH);

    // 16-bit or 32-bit endian correction
    if (s.endian) endianSwap(s, LENGTH, tmp0);
  
    char *tmpNoHead = new char [DATALENGTH];
    //get rid of headers and footers by memcpy 
    for (unsigned long j = 0; j<noFrames;j++) {
        memcpy(tmpNoHead + j*s.BYTESIZE_DET, tmp0 + s.HEADERSIZE + j*FRAMELENGTH, s.BYTESIZE_DET);
    }

  	//read whole line into a new 3D array
  	s.STEP = noFrames;
   	af::array ar_line_small = pointerToArray(s, tmpNoHead);
    delete[] tmp0;
    delete[] tmpNoHead;     
     
    //ar_line_small = ar_line_small * af::tile(generate_hann(ar_line_small),1,1,ar_line_small.dims(2));

    //we can do this as averaging trick (smear in x direction)
    //ar_line_small = ar_line_small + af::shift(ar_line_small,0,0,1)/2.0 + af::shift(ar_line_small,0,0,-1)/2.0;
    return ar_line_small;
}

// af::array get_disks_RANDOM(std::ifstream& fs, DefaultVals s,  unsigned long noFrames) {
//     //calculate size of read out data for the *tmp0 pointer in bytes
//     long int FRAMELENGTH = s.HEADERSIZE + s.BYTESIZE_DET + s.FOOTERSIZE;
//     long int LENGTH = noFrames * FRAMELENGTH;
//     long int DATALENGTH = s.BYTESIZE_DET * noFrames;
//     char *tmp0 = new char [LENGTH];
//     char *tmpNoHead = new char [DATALENGTH];
//   
//     //stack overflow answer to random number generation
//     std::mt19937 rng;
//     rng.seed(std::random_device()());
//     // distribution to search for random set of frames not over the edge of acquisition
//     unsigned long limit = s.SCAN_X;
//     if (limit>s.SCAN_Y) limit = s.SCAN_Y;
//     std::uniform_int_distribution<std::mt19937::result_type> dist6(0,limit - noFrames); 
//     
//     unsigned long rx = dist6(rng);
//     unsigned long ry = dist6(rng);
//     //somehow the random number overflows 
//     while (rx>limit-noFrames) rx = dist6(rng);
//     while (ry>limit-noFrames) ry = dist6(rng);
//     
//     unsigned long newpos = FRAMELENGTH * rx + FRAMELENGTH * ry * s.SCAN_X;
// 
//     fs.seekg(newpos, ios::beg);
// 
//     fs.read((char*)tmp0, LENGTH);
//     if (s.endian) endianSwap(s, LENGTH, tmp0);
// 
//     //get rid of headers and footers by memcpy 
//     for (int j = 0; j<noFrames;j++) {
//         memmove(tmpNoHead + j * s.BYTESIZE_DET , tmp0 + s.HEADERSIZE + j * FRAMELENGTH, s.BYTESIZE_DET );
//     }
//     
//     af::array ar_line = pointerToArray(s, tmpNoHead);
// 
//     delete[] tmp0;
//     delete[] tmpNoHead;  
//     
//     return af::resize(ar_line,int(ar_line.dims(0)*1.5),int(ar_line.dims(1)*1.5),AF_INTERP_BILINEAR);
// }

std::string sConvert (float number)
{
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

void endianSwap(DefaultVals s, int LENGTH, char* tmp) {
    //char*c0 = (char*)tmp0;
    //for (int i = 0; i < LENGTH*2; i += 2)
    //swap(c0[i], c0[i + 1]);
    // using GCC builtin swap because it is about 30% faster 
    if (s.bytedepth == 2) {
        unsigned short*c0 = (unsigned short*)tmp;
        for (int i = 0; i < LENGTH/2; i ++)
            c0[i] = __builtin_bswap16(c0[i]);
    }

    if (s.bytedepth == 4) {
        int*c0 = (int*)tmp;
        for (int i = 0; i < LENGTH/4; i ++)
            c0[i] = __builtin_bswap32(c0[i]);
    }
}

af::array pointerToArray(DefaultVals s, char* tmpNoHead, bool mask) {
    af::array ar_line_small;
    if (mask) s.STEP = 1;
    if (s.bytedepth == 1) {
        if (s.sign) ar_line_small = array(s.SIZE_DET,s.SIZE_DET,s.STEP, (char*) tmpNoHead);
        else ar_line_small= array(s.SIZE_DET,s.SIZE_DET,s.STEP, (unsigned char*) tmpNoHead);        
    }    
    if (s.bytedepth == 2) {
        if (s.sign) ar_line_small = array(s.SIZE_DET,s.SIZE_DET,s.STEP, (short*) tmpNoHead);
        else ar_line_small =  array(s.SIZE_DET,s.SIZE_DET,s.STEP, (unsigned short*) tmpNoHead);         
    }
    if (s.bytedepth == 4) {
        ar_line_small= array(s.SIZE_DET,s.SIZE_DET,s.STEP, (float*) tmpNoHead);
    }
    return ar_line_small;
}

af::array downscale(const DefaultVals &s,const array &ar_line_small) {
    return af::resize(ar_line_small,ar_line_small.dims(0)*s.DOWNSCALE,ar_line_small.dims(1)*s.DOWNSCALE,AF_INTERP_BILINEAR);
}



af::array zeropad(const DefaultVals &s, const af::array &ar_line_small) {
    unsigned long padding = ar_line_small.dims(0) / 4 * s.ZEROPAD;
    af::dim4 padn = af::dim4(padding,padding,0,0);
    return af::pad(ar_line_small,padn,padn,AF_PAD_ZERO);
}

//pre af3.7 zeropad functions
/*
af::array zeropad(const DefaultVals &s, const af::array &ar_line_small) {
    unsigned long newSize = ar_line_small.dims(0) * s.ZEROPAD;
    unsigned long L = (newSize - ar_line_small.dims(0)) / 2;
    unsigned long H = L + ar_line_small.dims(0) - 1 ;
    array ar_line(newSize,newSize,ar_line_small.dims(2),f32);
    ar_line=0.0;
    ar_line(seq(L,H),seq(L,H),seq())=ar_line_small.as(f32);
    return ar_line;
}
af::array zeropad(float scale, const af::array &ar_line_small) {
    unsigned long newSize = roundFloat(ar_line_small.dims(0) * scale); 
    unsigned long L = (newSize - ar_line_small.dims(0)) / 2;
    unsigned long H = L + ar_line_small.dims(0) - 1 ;
    array ar_line(newSize,newSize,ar_line_small.dims(2),f32);
    ar_line=0.0;
    ar_line(seq(L,H),seq(L,H),seq())=ar_line_small.as(f32);
    return ar_line;
}*/

af::array zeropad(float scale, const af::array &ar_line_small) {
    unsigned long padding = roundFloat(ar_line_small.dims(0) * scale / 2.0) - ar_line_small.dims(0) / 2 ;
    af::dim4 padn = af::dim4(padding,padding,0,0);
    return af::pad(ar_line_small,padn,padn,AF_PAD_ZERO);
}



array readDiskLine(ifstream& fs, DefaultVals s) {
    //calculate size of read out data for the *tmp0 pointer in bytes

    int FRAMELENGTH = s.HEADERSIZE+s.BYTESIZE_DET+s.FOOTERSIZE;
    int LENGTH = FRAMELENGTH * s.STEP; 
    int DATALENGTH = s.BYTESIZE_DET * s.STEP;

    char *tmp0 = new char [LENGTH];
    fs.read((char*)tmp0, LENGTH);

    // 16-bit or 32-bit endian correction
    if (s.endian) endianSwap(s, LENGTH, tmp0);
  
    char *tmpNoHead = new char [DATALENGTH];
    //get rid of headers and footers by memcpy 
    for (int j = 0; j<s.STEP;j++) {
        memcpy(tmpNoHead + j*s.BYTESIZE_DET, tmp0 + s.HEADERSIZE + j*FRAMELENGTH, s.BYTESIZE_DET);
    }

  	//read whole line into a new 3D array
   	af::array ar_line_small = pointerToArray(s, tmpNoHead);

    delete[] tmp0;
    delete[] tmpNoHead;     
     
    //ar_line_small = ar_line_small * af::tile(generate_hann(ar_line_small),1,1,ar_line_small.dims(2));

    //we can do this as averaging trick (smear in x direction)
    //ar_line_small = ar_line_small + af::shift(ar_line_small,0,0,1)/2.0 + af::shift(ar_line_small,0,0,-1)/2.0;
    //if (s.DOWNSCALE < 0.99) ar_line_small = downscale(s,ar_line_small);
    //if (s.ZEROPAD > 1) ar_line_small = zeropad(s,ar_line_small);
    return ar_line_small;
}

af::array generate_edge(af::array disk) {
     af::array gKerr,dx,dy;  
     gKerr = af::gaussianKernel(disk.dims(0), disk.dims(1),1,1);
     grad(dx,dy,gKerr);     
     array mask_dx = fftConvolve2(dx,disk,AF_CONV_DEFAULT );
     array mask_dy = fftConvolve2(dy,disk,AF_CONV_DEFAULT );
     array edge = af::hypot(mask_dx,mask_dy);
     return maxnorm(edge);
}

af::array generate_edge(DefaultVals s, af::array mask) {
    //edge creation
    if (s.DOG) {
        array gKerrL = af::gaussianKernel( mask.dims(0), mask.dims(1),s.dogLarge,s.dogLarge);
        array gKerrS = af::gaussianKernel( mask.dims(0), mask.dims(1),s.dogSmall,s.dogSmall);
        return af::fftConvolve2(gKerrL,mask)-af::fftConvolve2(gKerrS,mask);
    } else if (s.sigma > 0.19) {
        array gKerr,dx,dy;  
        gKerr = af::gaussianKernel( mask.dims(0), mask.dims(1),s.sigma,s.sigma);
        grad(dx,dy,gKerr);     

        array mask_dx = fftConvolve2(dx,mask,AF_CONV_DEFAULT );
        array mask_dy = fftConvolve2(dy,mask,AF_CONV_DEFAULT );
        return af::hypot(mask_dx,mask_dy);
    } else {
    // some datasets do really well with simple gradient (it probably requires smoothed mask done by AUTO algorithm)
    array dx,dy;
    grad(dx,dy,mask);     
    return af::hypot(dx,dy);
    }
}

// af::array generate_edge(DefaultVals s, af::array mask) {
//     //edge creation
//     if (s.DOG) {
//         array gKerrL = af::gaussianKernel( mask.dims(0), mask.dims(1),s.dogLarge,s.dogLarge);
//         array gKerrS = af::gaussianKernel( mask.dims(0), mask.dims(1),s.dogSmall,s.dogSmall);
//         return af::fftConvolve2(gKerrL,mask)-af::fftConvolve2(gKerrS,mask);
//     }
//     else{
//         
//         //cross correlation with upscale
//         array gKerr,dx,dy;  
//         gKerr = af::gaussianKernel( mask.dims(0), mask.dims(1),s.sigma,s.sigma);
//         grad(dx,dy,gKerr);     
// 
//         array MASK = af::fft2(mask,mask.dims(0)*2,mask.dims(1)*2);
//         array DX = af::fft2(dx,dx.dims(0)*2,dx.dims(1)*2);
//         array DY = af::fft2(dy,dy.dims(0)*2,dy.dims(1)*2);
//         array mask_dx = real(ifft2(MASK*DX));
//         array mask_dy = real(ifft2(MASK*DY));
//         
//         
//     
//         array edge = af::hypot(mask_dx,mask_dy);
// //edge = shift(edge,-edge.dims(0)/4,-edge.dims(1)/4);
// //    ccor_line /= float(ccor_line.dims(0))*float(ccor_line.dims(1));
//     unsigned int l = mask.dims(0)/2*(2 - 1);
//     unsigned int h = mask.dims(0)/2*(2 + 1) -1;
//     //return 
//     return edge(seq(l,h),seq(l,h));
//     }
// }

af::array generate_edge2(const af::array &disk, float smooth) {
     af::array gKerr,dx,dy;  
     unsigned int scale = 2;
     gKerr = af::gaussianKernel(disk.dims(0), disk.dims(1),smooth,smooth);
     grad(dx,dy,gKerr);     
     array mask_dx = fftConvolve2(dx,disk,AF_CONV_DEFAULT);
     array mask_dy = fftConvolve2(dy,disk,AF_CONV_DEFAULT);
     array edge = af::hypot(mask_dx,mask_dy);
     af::deviceGC();

     return edge;
}


af::array zoom_image(const af::array &disk, float zoom) {
    //create transformation matrix according to af documentation
    zoom = 1.0f/zoom; 
    float tx = float(disk.dims(0))*(1.0-zoom)*0.5;
    float ty = float(disk.dims(1))*(1.0-zoom)*0.5;
    float hA[] = {zoom, 0.0, tx, 0.0, zoom, ty};
    af::array trans(3,2,hA); 
    //zoom the array
    af::array zoomed = af::transform(disk,trans,disk.dims(0),disk.dims(1),AF_INTERP_BILINEAR);
    return zoomed;
}


void check_s(DefaultVals &s) {
    //dimension of computed arrays, number of parallel frames depends on this - ideally larger card will give result anyways

    s.SCAN_XY        = s.SCAN_X * s.SCAN_Y;
    s.bytedepth      = s.bitdepth / 8;
    cout << "bytedepth = "<<s.bytedepth<<endl;
    s.BYTESIZE_DET   = s.bytedepth * s.SIZE_DET * s.SIZE_DET;
}

// af::array get_mask(DefaultVals s) { 
//     //this is only programed for 16 bit mask - 8/32bits to be programmed
//     //READ TOP HAT MASK from maskfile
// 	ifstream th;
//     th.open(s.maskfile.c_str(), ios::in | ios::binary);
// 	
//     if (!th) {
// 	  cout << "\nSelected mask file does not exist!\n\n";
// 	  exit(1);
// 	}
//     char *tmp0 = new char [s.BYTESIZE_DET];
//     th.seekg(0);
//     th.read(tmp0, s.BYTESIZE_DET);
//     th.close();
// 
//     if (s.endian) endianSwap(s, s.BYTESIZE_DET, tmp0);
// 
//     //load mask pointer into mask array - small because it is not zeropadded
//     
//     array mask_small = pointerToArray(s,tmp0,true);
//     //create array for zeropadding and load mask into it
//     //if (s.DOWNSCALE < 0.99) mask_small = downscale(s,mask_small);
//     //if (s.ZEROPAD > 1) mask_small = zeropad(s,mask_small);
//     return mask_small.as(f32);
// }

void setDefaultSettings(DefaultVals &s) {
    //default values for analysis data needs to be 16bit unsigned - otherwise convert by imageJ or ..
    s.AUTO       = false;    // automatic edge/mask generation
    s.DOG        = false;	 // true for difference of gaussians to be used
    s.sigma      = 0.1;	   	 // sigma for dx,dy Gradient analysis
    s.STEP       = 256;		 // number of images analysed at same time (257 is prime number)
    s.HEADERSIZE = 384;	     // 384 new files and 256 old ones
    s.FOOTERSIZE = 0;	     // empad files have 1024 byte footers
    s.endian     = true;                 // correct endians - if used in windows or data were acquired on different endian machine
    s.show       = true;                 // show the disks - slow
    s.gradient   = false;
    s.dirname = "outputDir"; //dir_name and position of files
    
    //DEFAULT IMPORTANT DIMENSIONS OF SCAN ANALYSIS (adjust according the scan and needed 0 padding)
    s.ZEROPAD        = 2;     //doubles the size of the image by zeropadding to protect for close to edge effects in cross-correlation processing and Nyqist errors
    s.DOWNSCALE      = 1;    // for faster but less precise analysis
    s.SIZE_DET       = 256;     //256 for Medipix3 single
    s.SCAN_X         = 257;     //257
    s.SCAN_Y         = 256;	   //256

    s.sign           = false;    //switch for signed/unsigned data from detector
    s.bitdepth       = 16;
    //SIZE of the processing window!!
}


af::array partialConvolve(DefaultVals s, af::array edge_line,af::array EdgeLong) //done to speed up convolution of the edges - saving one fft which is done before sequential calculation
{   
    array Edge_line = fft2(edge_line,edge_line.dims(0)*s.ZEROPAD,edge_line.dims(1)*s.ZEROPAD);
    Edge_line = Edge_line*tile(EdgeLong,1,1,Edge_line.dims(2));
	array ccor_line = real(ifft2(Edge_line));
	// INVERT QUADRANTS AFTER IFFT
	ccor_line = shift(ccor_line,ccor_line.dims(0)/2,ccor_line.dims(1)/2);
    unsigned int l = edge_line.dims(0)/2*(s.ZEROPAD - 1);
    unsigned int h = edge_line.dims(1)/2*(s.ZEROPAD + 1) -1;
    return ccor_line(seq(l,h),seq(l,h));
}


af::array cross_correlate(const af::array &a1, const af::array &a2) 
{   
    array A1 = af::fft2(a1);
    array A2 = af::fft2(a2);
    array ccor_line = af::real(af::ifft2(A1*af::conjg(A2)));
	// INVERT QUADRANTS AFTER IFFT
	ccor_line = shift(ccor_line,ccor_line.dims(0)/2,ccor_line.dims(1)/2);

    return ccor_line/float(ccor_line.dims(0))/float(ccor_line.dims(1));
}

Results createResults(DefaultVals s) {
    Results calcs;
    calcs.sums = af::constant(0.0,s.SCAN_XY);
    calcs.rotsim = af::constant(0.0,s.SCAN_XY);
    calcs.rotmirror = af::constant(0.0,s.SCAN_XY);
    calcs.pacbed = af::constant(0.0,s.SIZE_DET,s.SIZE_DET );
    
    return calcs;
}
    
void generate_outputs(DefaultVals s,Results &calcs) {
    //sum and maximum of correlation image
    calcs.sums = af::moddims(calcs.sums, s.SCAN_X,s.SCAN_Y);
    
    calcs.rotsim = af::moddims(calcs.rotsim,s.SCAN_X,s.SCAN_Y);
    calcs.rotmirror = af::moddims(calcs.rotmirror,s.SCAN_X,s.SCAN_Y);
}

void data_check(ifstream& fs, DefaultVals s) {
	//check for datafile existence
	if (!fs) {
	    cout << "\nSelected raw/bin/mib/* file does not exist!\n\n"; exit(0);
	}
	
    fs.seekg(0, ios::end);

	ulong filesize = fs.tellg();
    
    cout << "header = " << s.HEADERSIZE << endl;
    cout << "footer = " << s.FOOTERSIZE << endl;
    cout << "detector = " << s.BYTESIZE_DET << endl;
    
    ulong dataframes = ulong(filesize/(s.HEADERSIZE+s.FOOTERSIZE+s.BYTESIZE_DET));
	cout << "\nSize of the file is: " << filesize << " B (" << filesize/1024/1024 << " MB / " << filesize/1024/1024/1024 <<"GB)";
	cout << "\nNumber of frames in binary file: " << dataframes;
	cout << "\nNumber of frames in xy parameters: " << s.SCAN_XY;
	  
    if (dataframes < s.SCAN_XY) {
	    cout << "\nFile is too small for the number of scan -x and -y parameters!\n\n";
	    exit(0);
	}
	
    //go back to begining
    fs.seekg(0, ios::beg);
}

int compute(ifstream& fs, DefaultVals &s, Results &calcs, af::Window& myWindow, bool end) {
    af::array ar_line, diskmask, z_array,z_rot_line,rot_line,ccor_line_rot,ccor_line_mirror;
    float deltaT = 1.0;
    
    while (s.frame_i < s.SCAN_XY - s.leftovers || end) {  
        end = false;
        
        //read sets of disk from the current position in fs stream
        ar_line = readDiskLine(fs,s); 

        //add PACBED contribution
        calcs.pacbed = calcs.pacbed + af::sum(ar_line,2);
        
        //generate low and high index for data storage
        long indL = s.frame_i;
        long indH = indL + s.STEP - 1;
        
        //sum data with a vector algorithm
        calcs.sums(seq(indL,indH)) = af::flat(af::sum(af::moddims(ar_line,ar_line.dims(0)*ar_line.dims(1),ar_line.dims(2)),0));
        ar_line = zeropad(s,ar_line);
        //remove noise hot pixels at the edge of gradient        
        diskmask = smoothed_disk(ar_line.dims(0),int(ar_line.dims(0)*0.45),5);
        ar_line = ar_line * 1.0f;
        ar_line = af::medfilt(ar_line);
        
        if (s.gradient) {
            ar_line = generate_edge(s,ar_line);
            ar_line = af::medfilt(ar_line);
        }
        
        ar_line = af::tile(diskmask,1,1,ar_line.dims(2)) * ar_line;
        normalizeCCinPlace(ar_line);

        myWindow.image(normaliseView(ar_line.slice(0)));
        
        z_array = zeropad(1.3,ar_line);
        unsigned long L = (z_array.dims(0) - ar_line.dims(0)) / 2;
        unsigned long H = L + ar_line.dims(0) - 1 ;
        
        //assumes square detector!!
        float angle = M_PI*float(s.angle)/180.0f;
        z_rot_line = af::rotate(z_array.copy(), angle,true, AF_INTERP_BILINEAR);
        
        // filter for random edge problem of interpolation
        rot_line = z_rot_line(seq(L,H),seq(L,H),seq());
        rot_line = af::tile(diskmask,1,1,rot_line.dims(2)) * rot_line;
                            
        normalizeCCinPlace(rot_line);
        
        ccor_line_rot = cross_correlate(ar_line,rot_line);  
        ccor_line_mirror = cross_correlate(rot_line,af::flip(rot_line,0)); 

        
        calcs.rotsim(seq(indL,indH)) = af::flat(af::max(af::max(ccor_line_rot,1),0));
        calcs.rotmirror(seq(indL,indH)) = af::flat(af::max(af::max(ccor_line_mirror,1),0));
    
        //deviceGC();
        //af::sync();
        
        s.frame_i+=s.STEP;
        
        // do some estimation for finishing time (hopefully better than windows LOL)
        float perc = 100.0*(float(s.frame_i)/float(s.SCAN_XY));
        float dt = timer::stop(s.start);
        float framerate = s.STEP /  (dt - deltaT);
        deltaT = dt;
        dt = dt / perc * 100.0 - dt;
        cout << "\r"<< roundFloat(perc) << "%      \t finishing in ~ "  << roundFloat(dt) << " s     \t" << roundFloat(framerate) << " fps      " << s.frame_i << " done out of " << s.SCAN_XY << "        ";

        //allow for termination by closing the window - if opened
        if (s.show and myWindow.close()) {
            cout << "\r=== CALCULATION TERMINATED === " << endl; 
            exit(1);
        }
    }
}


void debug_win(af::array a,af::array b, af::array c, af::array d,af::array e,af::array f,af::array g,af::array h, af::Window& myWindow) {
    myWindow.grid(2,4);
    myWindow(0,0).image(a, "sum");
    myWindow(0,1).image(b, "raw disk");
    myWindow(0,2).image(c, "edge");
    myWindow(0,3).image(d, "cross-corr");    
    myWindow(1,0).image(e, "dpcX");
    myWindow(1,1).image(f, "dpcY");
    myWindow(1,2).image(g, "corr-val");
    myWindow(1,3).image(h, "magnitude");

    myWindow.show();
}


void createLogFile(DefaultVals s) {
    //GENERATE LOGFILE RECORDING USED PARAMETERS
	string logfile="./log.txt";
    ofstream outputFile(logfile.c_str()); //this may not require c string conversion

    outputFile << "_____________________________________" << endl;
    outputFile << "__________ANALYSIS LOG FILE__________" << endl;
    outputFile << "_____________________________________" << endl;

    //outputFile  <<"mask file: "<< s.maskfile <<endl;
    outputFile  <<"binary file: "<< s.binfile <<endl;

    if (s.DOG) {
        outputFile  <<"difference of Gaussian filter used"<<endl;
        outputFile  <<"sigma1: "<< s.dogSmall <<endl;
        outputFile  <<"sigma2: "<< s.dogLarge <<endl;
    }
    else {
            outputFile  <<"Gaussian edge convolution filter used"<<endl;
            outputFile  <<"sigma: "<< s.sigma <<endl;
    }

    outputFile<<"output dir: "<< s.dirname <<endl; //this is fairly obvious
    outputFile<<"header size: "<< s.HEADERSIZE <<" kb" <<endl;
    outputFile<<"parallel: "<< s.STEP <<endl;
    outputFile<<"scan dimensions: X = "<< s.SCAN_X<<"  Y = "<< s.SCAN_Y <<endl;
    outputFile<<"zero padding: "<< s.ZEROPAD <<endl;
    outputFile<<"size of detector (square): "<< s.SIZE_DET <<endl;
    outputFile<<"calculation time (seconds): "<< timer::stop() <<endl;
    outputFile.close();
}

int getSettings(DefaultVals &s, int argc, char** argv) {
   while ((s.ind = getopt (argc, argv, "e:b:p:o:c:x:y:z:a:d:s:hwlfA:E")) != -1)
        switch (s.ind) {
            case 'h':
                showhelpinfo();
                return 1;
            case 'e':
                s.HEADERSIZE = atoi(optarg);
                break;
            case 'c':
                s.FOOTERSIZE = atoi(optarg);
                break;
            case 'p':
                s.STEP = atoi(optarg);
                break;
            case 'd':
                s.SIZE_DET = atoi(optarg);
                break;
            case 'x':
                s.SCAN_X = atoi(optarg);
                break;
            case 'y':
                s.SCAN_Y = atoi(optarg);
                break;
            case 'a':
                s.bitdepth = atoi(optarg);
                break;
            case 'z':
                s.ZEROPAD = atoi(optarg);
                break;
            case 'b':
                s.binfile = optarg;
                break;
            case 'o':
                s.dirname = optarg;
                break;        
            case 'w':
                s.endian = false;
                break;
            case 'A':
                s.angle = atof(optarg);
                break;            
            case 'f':
                s.show=false;
                break;
            case 'E':
                s.gradient=true;
                break;
            case 'l':
                s.sign=true; // if l than do unsigned calculation or some other binary shift
                break;
            case 's':
                s.sigma =  atof(optarg);
                break;
            default:
                showhelpinfo();
                return 1;
        }
}

void gotHere() {
    cout << "\nGot here!!" << endl;
    exit(0);
}

void saveimages(DefaultVals s, Results calcs) {
    //rotate the arrays - to correct detector to c++ scan rotation
    // shifts are changed due to the new pixel size if downscaled
    calcs.sums=rotate(calcs.sums,-M_PI/2.0f,false);
    calcs.rotsim=rotate(calcs.rotsim,-M_PI/2.0f,false);
    calcs.rotmirror=rotate(calcs.rotmirror,-M_PI/2.0f,false);

    //write the arrays

    //32 bit results
	saveImageNative("./pacbed.tif",calcs.pacbed);
	saveImageNative("./sum.tif",calcs.sums);
	saveImageNative("./rotation_symmetry.tif",calcs.rotsim);
	saveImageNative("./mirror_symmetry.tif",calcs.rotmirror);


    //16bit results
	calcs.sums=  normalizeImage(calcs.sums);
	calcs.pacbed = normalizeImage(calcs.pacbed);
    calcs.rotsim=  normalizeImage(calcs.rotsim);
	calcs.rotmirror = normalizeImage(calcs.rotmirror);

    saveImageNative("./16bit_sum.tif",calcs.sums.as(u16));
	saveImageNative("./16bit_pacbed.tif",calcs.pacbed.as(u16));
    saveImageNative("./16bit_rotation_symmetry.tif",calcs.rotsim.as(u16));
	saveImageNative("./16bit_mirror_symmetry.tif",calcs.rotmirror.as(u16));

    createLogFile(s);
}


void showMe() {
    
    printf(R"EOF(
 __                     ______    
(_   ._ _ ._ _  __|_._ (_  ||_|\/|
__)\/| | || | |(/_|_|\/__) ||_|  |
   /                 /            

)EOF");
	cout << "\n\n (C) Matus Krajnak\n\n";
}

void calcString() {
    cout << "\n\n=========================\n";
    cout << "== CALCULATION STARTED ==\n";
    cout << "=========================\n\n";
}
