#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <stdint.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "algo_top.h"
#include "algo_top_parameters.h"

using namespace std;

int main() {

	ap_uint<576> link_in[N_INPUT_LINKS];
	ap_uint<576> link_out[N_OUTPUT_LINKS];

    size_t start = 0;
    size_t end = 13;

    for(size_t i=0; i<N_INPUT_LINKS; i++){
    	   link_in[i] = 0;
       }
    //sector 2
    link_in[0]  = "0x00000000000000000000000000000000000000000000000000000000000E00000000000000000000000000000000000000000000000000000000000000000000000000000000000A";
    link_in[1]  = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[2]  = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[3]  = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[4]  = "0x000000000000000000000000000000000000000000000000000000000000007800000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[5]  = "0x00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000B";
    link_in[6]  = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[7]  = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[8]  = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[9]  = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[10] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[11] = "0x000000000000000000000000000000000000000000000000000000000000000000000000050000000000000000000000000000000000000000000000003000000000000000000000";
    link_in[12] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[13] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[14] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[15] = "0x00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000003C000";

    //sector 4
    link_in[16] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[17] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[18] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[19] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[20] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000240000000000000000000000000000";
    link_in[21] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[22] = "0x0000000000000000000000000000000000000000000000000000000000000000000002C0000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[23] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[24] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[25] = "0x00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000018000000000000000000000000000000C";
    link_in[26] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[27] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[28] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[29] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[30] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    link_in[31] = "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";

	// Run the algorithm

	  algo_top(link_in, link_out);

	  	cout << hex << "link_out[0]: " << link_out[0] << endl;
	  	cout << hex << "link_out[1]: " << link_out[1] << endl;
//	  	cout << hex << "link_out[2]: " << link_out[2] << endl;
//	  	cout << hex << "link_out[3]: " << link_out[3] << endl;

	  cout << "printing towers et" << endl;

	  		start = 0;
	  		for(loop oLink=0; oLink<32; oLink++){
	  			size_t end = start + 11;
	  			cout << link_out[0].range(end, start) << " " ;
                                if(oLink%4 == 3) cout << endl ;
	  			start += 12;
	  		}

	  		cout << "printing cluster et" << endl;

for(loop i=0; i<301; i=i+60)
cout << link_out[1].range(11+i, 0+i) << " " << link_out[1].range(16+i, 12+i) << " " << link_out[1].range(18+i, 17+i) << " " << link_out[1].range(21+i, 19+i) << " " << link_out[1].range(24+i, 22+i) << " " << link_out[1].range(27+i, 25+i) << " " << link_out[1].range(42+i, 28+i) << " " << link_out[1].range(57+i, 43+i) << " " << link_out[1].range(59+i, 58+i) <<  endl ;

return 0;

}



