#include "symmetry_support.h"

using namespace std;
using namespace af;


//main function
int main(int argc, char** argv ) {
    
    // who made this
    showMe();
    // tell what background arrayfire is going to use
    cout << "\n"; af::info();

    // if no options given show help and exit
    if(argc == 1) {
        showhelpinfo();
        exit(1);
    }
    
    //create and check a struct for setup of s - everything can be passed to functions ;) all defined in dpc_support.h/.cpp
    DefaultVals s; 
    setDefaultSettings(s);    
    // set the important variables by passing onto function from dpc_support.cpp and check the basics for running the analysis
    getSettings(s,argc,argv);
    check_s(s);

    //OPEN BINARY BIG FILE FOR ANALYSIS - open at end for filelength check
	ifstream fs;
    fs.open(s.binfile.c_str(), ios::out | ios::binary | ios::ate);
        
    //use hanning window - its probably wrongly implemented
    s.hann = false;
    
    //file checks and go back to beginning
    data_check(fs, s);
    
    Window myWindow(256*4, 64 + 256*2, "computation visualisation");
    if (!s.show) myWindow.setVisibility(false);

    //read data by for loop one s.STEP at a time
    s.leftovers = s.SCAN_XY % s.STEP;

    //show that calculation has started
    calcString();
    //initiate Results
    Results calcs;
    //generate check for exception 
    bool done = false;

    //time operational part of the code
    s.start = timer::start();

    //this while / try / catch is done to correct for running out of memory on gpu
    while (!done) {
        try {
            //go to the beginning of the file and restart frame counter
            fs.seekg(0, ios::beg);
            s.frame_i = 0;
            //arrays to store resulting calculations
            calcs = createResults(s);
            
            
            // ======CALCULATION START=======
            //do the bulk calculation
            compute(fs, s, calcs, myWindow, false);
            //deal with any leftovers at the end of the file if XY % STEP > 0
            if (s.leftovers > 0) {
                s.STEP = s.leftovers;
                compute(fs, s, calcs, myWindow, true);
            }
            // ======CALCULATION FINISH=======
            done = true;
            
        } catch (af::exception& ae) {
            if (ae.err() == 101 || 998) {
                s.STEP=int(s.STEP/2); // it is actually faster than / 2 
                s.leftovers = s.SCAN_XY % s.STEP;
                af::deviceGC();
                cout << "\nAF error 101/998: device run out of memory\n";
                cout << "trying with parralel frames halved to: " << s.STEP << endl;
            } else {
                std::cerr << ae.what() << std::endl;
                done = true;
             }
        }
    }
    //gotHere();
    float finito = timer::stop();
    cout << endl;
    printf("\relapsed seconds: %g \n", roundFloat(finito));
    cout << "average fps " << roundFloat(s.SCAN_XY / finito) << endl<<endl;
    
    //check if the directory exists and if yes add 0 at the end to avoid the loss of data. When created, chdir into it
    while (DirectoryExists(s.dirname.c_str())) s.dirname+="0";
	mkdir(s.dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    chdir(s.dirname.c_str());
    
    //BUILD RESULTING IMAGES
    generate_outputs(s,calcs);

    //save and show resulting images
    saveimages(s,calcs);

    fs.close();
    
    return 0;
}
