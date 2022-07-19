#include "algo_top.h"
#include "algo_top_parameters.h"
#include "bitonicSort16.h"

void processInputLinks(ap_uint<576> link_in[N_INPUT_LINKS], crystal ECALRegion3x4_1[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], crystal ECALRegion3x4_2[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI],crystal ECALRegion3x4_3[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI]){

                #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
        #pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_1 complete dim=0
        #pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_2 complete dim=0
        #pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_3 complete dim=0

        ap_uint<32> start = 0;
        ap_uint<32> end = 13;

        ap_uint<6> wordId, wordId1, wordId2, startId ;


        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
                //#pragma HLS UNROLL
                for(loop j=0; j<CRYSTAL_IN_PHI; j++){
                        #pragma HLS UNROLL
                       wordId = (i/5)*4+(j/5) ;
                       wordId1 = wordId + 12 ;
                       wordId2 = wordId + 24 ;
                       startId = (i%5)*5+(j%5) ;
                       start = startId*14 ; end = start + 13 ;

                       ECALRegion3x4_1[i][j] = crystal(link_in[wordId].range(end, start));
                       ECALRegion3x4_2[i][j] = crystal(link_in[wordId1].range(end, start));
                       if(wordId2 <= 31) ECALRegion3x4_3[i][j] = crystal(link_in[wordId2].range(end, start));
                       else ECALRegion3x4_3[i][j] = crystal(0);
                }
        }


}

void processOutLink(Cluster inCluster[3], tower_t inTower[16], ap_uint<576> &outLink){
        #pragma HLS PIPELINE II=9
        #pragma HLS LATENCY min=9

        ap_uint<32> start = 0;

        for(loop oLink=0; oLink<16; oLink++){
                #pragma HLS UNROLL
                ap_uint<32> end = start + 11;
                outLink.range(end, start) = inTower[oLink].et();
                start += 12;
        }

        for(loop oLink=0; oLink<3; oLink++){
                #pragma HLS UNROLL
                ap_uint<32> end = start + 59;
                outLink.range(end, start) = inCluster[oLink];
                start += 60;
        }

}

ecaltp_t bestOf2(const ecaltp_t& ecaltp0, const ecaltp_t& ecaltp1) {
        ecaltp_t x;
        x = (ecaltp0.energy > ecaltp1.energy)?ecaltp0:ecaltp1;

        return x;
}

ecaltp_t getPeakBin20N(const etaStrip_t& etaStrip){
//#pragma HLS PIPELINE II=2
#pragma HLS latency min=4

ecaltp_t best01 = bestOf2(etaStrip.cr0,etaStrip.cr1) ;
ecaltp_t best23 = bestOf2(etaStrip.cr2,etaStrip.cr3) ;
ecaltp_t best45 = bestOf2(etaStrip.cr4,etaStrip.cr5) ;
ecaltp_t best67 = bestOf2(etaStrip.cr6,etaStrip.cr7) ;
ecaltp_t best89 = bestOf2(etaStrip.cr8,etaStrip.cr9) ;
ecaltp_t best1011 = bestOf2(etaStrip.cr10,etaStrip.cr11) ;
ecaltp_t best1213 = bestOf2(etaStrip.cr12,etaStrip.cr13) ;
ecaltp_t best1415 = bestOf2(etaStrip.cr14,etaStrip.cr15) ;
ecaltp_t best1617 = bestOf2(etaStrip.cr16,etaStrip.cr17) ;
ecaltp_t best1819 = bestOf2(etaStrip.cr18,etaStrip.cr19) ;

ecaltp_t best0123 = bestOf2(best01,best23) ;
ecaltp_t best4567 = bestOf2(best45,best67) ;
ecaltp_t best891011 = bestOf2(best89,best1011) ;
ecaltp_t best12131415 = bestOf2(best1213,best1415) ;
ecaltp_t best16171819 = bestOf2(best1617,best1819) ;

ecaltp_t best01234567 = bestOf2(best0123,best4567) ;
ecaltp_t best89101112131415 = bestOf2(best891011,best12131415) ;

ecaltp_t best0to15 = bestOf2(best01234567,best89101112131415) ;
ecaltp_t bestOf20 = bestOf2(best0to15,best16171819) ;

return bestOf20 ;
}

crystalMax getPeakBin15N(const etaStripPeak_t& etaStrip){
//#pragma HLS PIPELINE II=2
#pragma HLS latency min=4

crystalMax x;

ecaltp_t best01 = bestOf2(etaStrip.pk0,etaStrip.pk1) ;
ecaltp_t best23 = bestOf2(etaStrip.pk2,etaStrip.pk3) ;
ecaltp_t best45 = bestOf2(etaStrip.pk4,etaStrip.pk5) ;
ecaltp_t best67 = bestOf2(etaStrip.pk6,etaStrip.pk7) ;
ecaltp_t best89 = bestOf2(etaStrip.pk8,etaStrip.pk9) ;
ecaltp_t best1011 = bestOf2(etaStrip.pk10,etaStrip.pk11) ;
ecaltp_t best1213 = bestOf2(etaStrip.pk12,etaStrip.pk13) ;

ecaltp_t best0123 = bestOf2(best01,best23) ;
ecaltp_t best4567 = bestOf2(best45,best67) ;
ecaltp_t best891011 = bestOf2(best89,best1011) ;
ecaltp_t best121314 = bestOf2(best1213,etaStrip.pk14) ;

ecaltp_t best01234567 = bestOf2(best0123,best4567);
ecaltp_t best891011121314 = bestOf2(best891011,best121314) ;

ecaltp_t bestOf15 = bestOf2(best01234567,best891011121314) ;

        x.energy = bestOf15.energy ;
        x.etaMax = bestOf15.eta ;
        x.phiMax = bestOf15.phi ;

return x ;
}

void getTowerEt(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<12> towerEt[12]){
//        #pragma HLS PIPELINE II=9
        #pragma HLS latency min=4
        #pragma HLS ARRAY_PARTITION variable=tempX complete dim=0
        #pragma HLS ARRAY_PARTITION variable=towerEt complete dim=0

        ap_uint<10>  temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI] ;
        #pragma HLS ARRAY_PARTITION variable=temp complete dim=0

        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
                      #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
            temp[i][k] = tempX[i][k].energy ;
         }}

         for(loop i=0; i<CRYSTAL_IN_ETA; i=i+5){
              #pragma HLS UNROLL
              for(loop k=0; k<CRYSTAL_IN_PHI; k=k+5){
                 #pragma HLS UNROLL
                 ap_uint<12> cntr = (i/5)*4+(k/5) , a, b, c, d, e;
                 a = temp[i][k] + temp[i][k+1] + temp[i][k+2] + temp[i][k+3] + temp[i][k+4] ;
                 b = temp[i+1][k] + temp[i+1][k+1] + temp[i+1][k+2] + temp[i+1][k+3] + temp[i+1][k+4] ;
                 c = temp[i+2][k] + temp[i+2][k+1] + temp[i+2][k+2] + temp[i+2][k+3] + temp[i+2][k+4] ;
                 d = temp[i+3][k] + temp[i+3][k+1] + temp[i+3][k+2] + temp[i+3][k+3] + temp[i+3][k+4] ;
                 e = temp[i+4][k] + temp[i+4][k+1] + temp[i+4][k+2] + temp[i+4][k+3] + temp[i+4][k+4] ;
                 towerEt[cntr]= a + b + c + d + e ;
          }}

}

clusterInfo getClusterPosition(const ecalRegion_t& ecalRegion){
 #pragma HLS latency min=4

        etaStripPeak_t etaStripPeak;
        clusterInfo cluster ;


                etaStripPeak.pk0 = getPeakBin20N(ecalRegion.etaStrip0);
                etaStripPeak.pk1 = getPeakBin20N(ecalRegion.etaStrip1);
                etaStripPeak.pk2 = getPeakBin20N(ecalRegion.etaStrip2);
                etaStripPeak.pk3 = getPeakBin20N(ecalRegion.etaStrip3);
                etaStripPeak.pk4 = getPeakBin20N(ecalRegion.etaStrip4);
                etaStripPeak.pk5 = getPeakBin20N(ecalRegion.etaStrip5);
                etaStripPeak.pk6 = getPeakBin20N(ecalRegion.etaStrip6);
                etaStripPeak.pk7 = getPeakBin20N(ecalRegion.etaStrip7);
                etaStripPeak.pk8 = getPeakBin20N(ecalRegion.etaStrip8);
                etaStripPeak.pk9 = getPeakBin20N(ecalRegion.etaStrip9);
                etaStripPeak.pk10 = getPeakBin20N(ecalRegion.etaStrip10);
                etaStripPeak.pk11 = getPeakBin20N(ecalRegion.etaStrip11);
                etaStripPeak.pk12 = getPeakBin20N(ecalRegion.etaStrip12);
                etaStripPeak.pk13 = getPeakBin20N(ecalRegion.etaStrip13);
                etaStripPeak.pk14 = getPeakBin20N(ecalRegion.etaStrip14);

        crystalMax peakIn15;
        peakIn15 = getPeakBin15N(etaStripPeak);

        cluster.seedEnergy = peakIn15.energy ;
        cluster.energy = 0 ;
        cluster.etaMax = peakIn15.etaMax ;
        cluster.phiMax = peakIn15.phiMax ;
        cluster.brems = 0 ;
        cluster.et5x5 = 0 ;
        cluster.et2x5 = 0 ;

return cluster ;
}

Cluster packCluster(ap_uint<15>& clusterEt, ap_uint<5>& etaMax_t, ap_uint<5>& phiMax_t, ap_uint<15>& et5x5_t, ap_uint<15>& et2x5_t, ap_uint<2>& brems_t ){

        ap_uint<12> peggedEt;
        Cluster pack;

        ap_uint<5> towerEta = (etaMax_t)/5;
        ap_uint<2> towerPhi = (phiMax_t)/5;
        ap_uint<3> clusterEta = etaMax_t - 5*towerEta;
        ap_uint<3> clusterPhi = phiMax_t - 5*towerPhi;
        ap_uint<15> clusterEt5x5 = et5x5_t ;
        ap_uint<15> clusterEt2x5 = et2x5_t ;
        ap_uint<2> clusterBrems = brems_t ;

        peggedEt = (clusterEt > 0xFFF)? (ap_uint<12>)0xFFF : (ap_uint<12>) clusterEt;

        pack = Cluster(peggedEt, towerEta, towerPhi, clusterEta, clusterPhi, 0, clusterEt5x5, clusterEt2x5, clusterBrems);

return pack;
}

void RemoveTmp(crystal temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi, ap_uint<2> brems  ){
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0
#pragma HLS latency min=3
//#pragma HLS PIPELINE II=9

        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
//                        #pragma HLS UNROLL
           if(i>=seed_eta-1 && i<=seed_eta+1){
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
            if(k>=seed_phi-2 && k<=seed_phi+2)  temp[i][k].energy = 0 ;}
         }}

        if(brems == 1){
        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
//                        #pragma HLS UNROLL
           if(i>=seed_eta-1 && i<=seed_eta+1){
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
            if(k>=seed_phi-2-5 && k<=seed_phi+2-5)  temp[i][k].energy = 0 ;}
                        }
                }
         }

        if(brems == 2){
        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
//                        #pragma HLS UNROLL
           if(i>=seed_eta-1 && i<=seed_eta+1 ){
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
            if(k>=seed_phi-2+5 && k<=seed_phi+2+5)  temp[i][k].energy = 0 ;}
         }}}
}

clusterInfo getBremsValuesPos(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi ){
#pragma HLS ARRAY_PARTITION variable=tempX complete dim=0
#pragma HLS latency min=6
//#pragma HLS PIPELINE II=9

        ap_uint<12> temp[CRYSTAL_IN_ETA+2][CRYSTAL_IN_PHI+4] ;
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0

        ap_uint<12> eta_slice[3] ;
#pragma HLS ARRAY_PARTITION variable=eta_slice complete dim=0

clusterInfo cluster_tmp;

        for(loop i=0; i<CRYSTAL_IN_ETA+2; i++){
                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI+4; k++){
                        #pragma HLS UNROLL
            temp[i][k] = 0 ;
         }}

        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
                        #pragma HLS UNROLL
           for(loop k=5; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
            temp[i+1][k-3] = tempX[i][k].energy ;
         }}


        ap_uint<6> seed_eta1,  seed_phi1 ;

        seed_eta1 = seed_eta ; //to start from corner
        seed_phi1 = seed_phi ; //to start from corner
// now we are in the left bottom corner
        ap_uint<12> tmp1, tmp2, tmp3 ;

        for(loop j=0; j<CRYSTAL_IN_ETA; j++){
//                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
              if(j== seed_eta1 && k == seed_phi1)
                 {
                for(loop m=0; m<3 ; m++){
                        #pragma HLS UNROLL
                tmp1 = temp[j+m][k] + temp[j+m][k+1] ;
                tmp2 = temp[j+m][k+2] + temp[j+m][k+3] ;
                tmp3 = tmp1 + temp[j+m][k+4] ;
                eta_slice[m] = tmp2 + tmp3 ;
//                eta_slice[m] = temp[j+m][k] + temp[j+m][k+1] +temp[j+m][k+2] +temp[j+m][k+3] +temp[j+m][k+4] ;
                        }
               }
          }}

         cluster_tmp.energy=eta_slice[0] + eta_slice[1] + eta_slice[2] ;

return cluster_tmp ;
}

clusterInfo getBremsValuesNeg(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi ){
#pragma HLS ARRAY_PARTITION variable=tempX complete dim=0
#pragma HLS latency min=6
//#pragma HLS PIPELINE II=9

        ap_uint<12> temp[CRYSTAL_IN_ETA+2][CRYSTAL_IN_PHI+4] ;
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0

        ap_uint<12> eta_slice[3] ;
#pragma HLS ARRAY_PARTITION variable=eta_slice complete dim=0

clusterInfo cluster_tmp;

        for(loop i=0; i<CRYSTAL_IN_ETA+2; i++){
                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI+4; k++){
                        #pragma HLS UNROLL
            temp[i][k] = 0 ;
         }}

        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI-5; k++){
                        #pragma HLS UNROLL
            temp[i+1][k+7] = tempX[i][k].energy ;
         }}


        ap_uint<6> seed_eta1,  seed_phi1 ;

        seed_eta1 = seed_eta ; //to start from corner
        seed_phi1 = seed_phi ; //to start from corner
// now we are in the left bottom corner

        ap_uint<12> tmp1, tmp2, tmp3 ;

        for(loop j=0; j<CRYSTAL_IN_ETA; j++){
//                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
              if(j== seed_eta1 && k == seed_phi1)
                 {
                for(loop m=0; m<3 ; m++){
                        #pragma HLS UNROLL
                tmp1 = temp[j+m][k] + temp[j+m][k+1] ;
                tmp2 = temp[j+m][k+2] + temp[j+m][k+3] ;
                tmp3 = tmp1 + temp[j+m][k+4] ;
                eta_slice[m] = tmp2 + tmp3 ;
                        }
               }
          }}

         cluster_tmp.energy=eta_slice[0] + eta_slice[1] + eta_slice[2] ;

return cluster_tmp ;
}


clusterInfo getClusterValues(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi ){
#pragma HLS ARRAY_PARTITION variable=tempX complete dim=0
#pragma HLS latency min=6
//#pragma HLS PIPELINE II=9

        ap_uint<12> temp[CRYSTAL_IN_ETA+4][CRYSTAL_IN_PHI+4] ;
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0

        ap_uint<12> eta_slice[5] ;
#pragma HLS ARRAY_PARTITION variable=eta_slice complete dim=0


        ap_uint<12> et2x5_1Tot, et2x5_2Tot, etSum2x5 ;
        ap_uint<12> et5x5Tot ;

clusterInfo cluster_tmp;

        for(loop i=0; i<CRYSTAL_IN_ETA+4; i++){
                       #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI+4; k++){
                        #pragma HLS UNROLL
            temp[i][k] = 0 ;
         }}

        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
            temp[i+2][k+2] = tempX[i][k].energy ;
         }}


        ap_uint<6> seed_eta1,  seed_phi1 ;

        seed_eta1 = seed_eta ; //to start from corner
        seed_phi1 = seed_phi ; //to start from corner
// now we are in the left bottom corner
        ap_uint<12> tmp1, tmp2, tmp3 ;

        for(loop j=0; j<CRYSTAL_IN_ETA; j++){
//                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
              if(j== seed_eta1 && k == seed_phi1)
                 {
                for(loop m=0; m<5 ; m++){
                        #pragma HLS UNROLL
                tmp1 = temp[j+m][k] + temp[j+m][k+1] ;
                tmp2 = temp[j+m][k+2] + temp[j+m][k+3] ;
                tmp3 = tmp1 + temp[j+m][k+4] ;
                eta_slice[m] = tmp2 + tmp3 ;
//                eta_slice[m] = temp[j+m][k] + temp[j+m][k+1] +temp[j+m][k+2] +temp[j+m][k+3] +temp[j+m][k+4] ;
                        }
               }
          }}


         cluster_tmp.energy=eta_slice[1] + eta_slice[2] + eta_slice[3] ;

          et5x5Tot = eta_slice[0] + eta_slice[1] + eta_slice[2] + eta_slice[3] + eta_slice[4] ;
          et2x5_1Tot = eta_slice[1] + eta_slice[2] ;
          et2x5_2Tot = eta_slice[2] + eta_slice[3] ;


          if(et2x5_1Tot >= et2x5_2Tot) etSum2x5 = et2x5_1Tot ;
          else etSum2x5 = et2x5_2Tot ;

          cluster_tmp.et5x5 = et5x5Tot ;
          cluster_tmp.et2x5 = etSum2x5 ;


return cluster_tmp ;
}


Cluster getRegion3x4(crystal temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> eta_offset){
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0
//#pragma HLS PIPELINE II=2
#pragma HLS latency min=24

Cluster returnCluster, forCluster;
clusterInfo cluster_tmp;
clusterInfo cluster_tmpCenter;
clusterInfo cluster_tmpBneg;
clusterInfo cluster_tmpBpos;

ecalRegion_t ecalRegion;

        ecalRegion = initStructure(temp) ;

        cluster_tmp = getClusterPosition(ecalRegion) ;

        ap_uint<5> seed_phi = cluster_tmp.phiMax ;
        ap_uint<5> seed_eta = cluster_tmp.etaMax ;

        cluster_tmpCenter = getClusterValues(temp, seed_eta, seed_phi) ;
        cluster_tmpBneg = getBremsValuesNeg(temp, seed_eta, seed_phi) ;
        cluster_tmpBpos = getBremsValuesPos(temp, seed_eta, seed_phi) ;

        cluster_tmp.energy = cluster_tmpCenter.energy;

        cluster_tmp.brems = 0 ;

         if(cluster_tmpBneg.energy > cluster_tmpCenter.energy/8 && cluster_tmpBneg.energy > cluster_tmpBpos.energy) {
            cluster_tmp.energy = cluster_tmpCenter.energy + cluster_tmpBneg.energy ; cluster_tmp.brems = 1 ; }
         else if(cluster_tmpBpos.energy > cluster_tmpCenter.energy/8){
            cluster_tmp.energy = cluster_tmpCenter.energy + cluster_tmpBpos.energy ; cluster_tmp.brems = 2 ; }

//eta, phi, seed energy in cluster_tmp; energy and brems in cluster_tmp1

        forCluster = packCluster(cluster_tmp.energy, cluster_tmp.etaMax, cluster_tmp.phiMax, cluster_tmp.et5x5, cluster_tmp.et2x5, cluster_tmp.brems);

        RemoveTmp(temp, seed_eta, seed_phi, cluster_tmp.brems ) ;

        ap_uint<5> towerEta = forCluster.towerEta() + eta_offset ;
        returnCluster = Cluster(forCluster.clusterEnergy(), towerEta, forCluster.towerPhi(), forCluster.clusterEta(), forCluster.clusterPhi(), forCluster.satur(), forCluster.et5x5(), forCluster.et2x5(), forCluster.brems());

return returnCluster;
}

void stitchClusters(Cluster clusterIn[Nbclusters], Cluster clusterOut[Nbclusters]){
        #pragma HLS ARRAY_PARTITION variable=clusterIn complete dim=0
        #pragma HLS ARRAY_PARTITION variable=clusterOut complete dim=0
        #pragma HLS latency min=10
//        #pragma HLS PIPELINE II=9
// we stitch 3 regions - 2 boundaries, the middle 3x4 region eta:3-5/12-14
// has lower 2:3/11:12 and upper 5:6/14:15 boundaries, we check lower and upper
// 3x4 regions and stitch if needed. The low et cluster gets et=0
//
       for(int i=0; i<Nbclusters; i++){
        #pragma HLS UNROLL
         clusterOut[i] = clusterIn[i] ;
         if(i >= 12) clusterOut[i] = Cluster(0,0,0,0,0,0,0,0,0) ;
       }

       for(int i=4; i<8; i++){
//      #pragma HLS UNROLL
        if( clusterIn[i].clusterEnergy() > 0 && clusterIn[i].towerEta() == 12 && clusterIn[i].clusterEta() == 0){
         for(int k=0; k<4; k++){
        #pragma HLS UNROLL
          if(clusterIn[k].clusterEnergy() > 0 && clusterIn[k].towerEta() == 11 && clusterIn[k].clusterEta() == 4){
          ap_uint<5> phi1 = clusterIn[i].towerPhi()*5 + clusterIn[i].clusterPhi() ;
          ap_uint<5> phi2 = clusterIn[k].towerPhi()*5 + clusterIn[k].clusterPhi() ;
          ap_uint<5> dPhi ; dPhi=(phi1 > phi2)?(phi1-phi2):(phi2-phi1) ;
           if( dPhi < 2 ) {
           ap_uint<12> one = clusterIn[i].clusterEnergy() ;
           ap_uint<12> two = clusterIn[k].clusterEnergy() ;
           ap_uint<12> sum = one+two ;
           ap_uint<15> et5x5 = clusterIn[i].et5x5() + clusterIn[k].et5x5() ;
           ap_uint<15> et2x5 = clusterIn[i].et2x5() + clusterIn[k].et2x5() ;
              if (one > two){
 clusterOut[i] = Cluster(sum, clusterIn[i].towerEta(), clusterIn[i].towerPhi(), clusterIn[i].clusterEta(), clusterIn[i].clusterPhi(), clusterIn[i].satur(), et5x5, et2x5, clusterIn[i].brems()) ;
 clusterOut[k] = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0) ;
                            }
              else {
 clusterOut[k] = Cluster(sum, clusterIn[k].towerEta(), clusterIn[k].towerPhi(), clusterIn[k].clusterEta(), clusterIn[k].clusterPhi(), clusterIn[k].satur(), et5x5, et2x5, clusterIn[k].brems()) ;
 clusterOut[i] = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0) ;
                   }
                 }}
          }
         }

        if( clusterIn[i].clusterEnergy() > 0 && clusterIn[i].towerEta() == 14 && clusterIn[i].clusterEta() == 4){
         for(int k=8; k<12; k++){
        #pragma HLS UNROLL
          if(clusterIn[k].clusterEnergy() > 0 && clusterIn[k].towerEta() == 15 && clusterIn[k].clusterEta() == 0){
          ap_uint<5> phi1 = clusterIn[i].towerPhi()*5 + clusterIn[i].clusterPhi() ;
          ap_uint<5> phi2 = clusterIn[k].towerPhi()*5 + clusterIn[k].clusterPhi() ;
          ap_uint<5> dPhi ; dPhi=(phi1 > phi2)?(phi1-phi2):(phi2-phi1) ;
           if( dPhi < 2 ) {
           ap_uint<12> one = clusterIn[i].clusterEnergy() ;
           ap_uint<12> two = clusterIn[k].clusterEnergy() ;
           ap_uint<12> sum = one+two ;
           ap_uint<15> et5x5 = clusterIn[i].et5x5() + clusterIn[k].et5x5() ;
           ap_uint<15> et2x5 = clusterIn[i].et2x5() + clusterIn[k].et2x5() ;
              if (one > two){
 clusterOut[i] = Cluster(sum, clusterIn[i].towerEta(), clusterIn[i].towerPhi(), clusterIn[i].clusterEta(), clusterIn[i].clusterPhi(), clusterIn[i].satur(), et5x5, et2x5, clusterIn[i].brems()) ;
 clusterOut[k] = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0) ;
                            }
              else {
 clusterOut[k] = Cluster(sum, clusterIn[k].towerEta(), clusterIn[k].towerPhi(), clusterIn[k].clusterEta(), clusterIn[k].clusterPhi(), clusterIn[k].satur(), et5x5, et2x5, clusterIn[k].brems()) ;
 clusterOut[i] = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0) ;
                   }
                 }}
          }}
        }

}

void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]){
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=9
#pragma HLS INTERFACE ap_ctrl_hs port=return

#pragma HLS latency min=100

crystal ECALRegion3x4_1[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_1 complete dim=0

crystal ECALRegion3x4_2[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_2 complete dim=0

crystal ECALRegion3x4_3[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_3 complete dim=0

//creating 2 15x20 crystals temporary and temporary1

processInputLinks(link_in, ECALRegion3x4_1, ECALRegion3x4_2, ECALRegion3x4_3);

Cluster sort_clusterIn[Nbclusters];
Cluster sort_clusterInStitched[Nbclusters];
Cluster sort_clusterOut[Nbclusters];

tower_t towerEt[36];

#pragma HLS ARRAY_PARTITION variable=sort_clusterIn complete dim=0
#pragma HLS ARRAY_PARTITION variable=sort_clusterInStitched complete dim=0
#pragma HLS ARRAY_PARTITION variable=sort_clusterOut complete dim=0
#pragma HLS ARRAY_PARTITION variable=towerEt complete dim=0

// 9 is eta offset, since this is second part of RCT card
//
        sort_clusterIn[0] = getRegion3x4(ECALRegion3x4_1, 9);
        sort_clusterIn[1] = getRegion3x4(ECALRegion3x4_1, 9);
        sort_clusterIn[2] = getRegion3x4(ECALRegion3x4_1, 9);
        sort_clusterIn[3] = getRegion3x4(ECALRegion3x4_1, 9);

        sort_clusterIn[4] = getRegion3x4(ECALRegion3x4_2, 12);
        sort_clusterIn[5] = getRegion3x4(ECALRegion3x4_2, 12);
        sort_clusterIn[6] = getRegion3x4(ECALRegion3x4_2, 12);
        sort_clusterIn[7] = getRegion3x4(ECALRegion3x4_2, 12);

        sort_clusterIn[8] = getRegion3x4(ECALRegion3x4_3, 15);
        sort_clusterIn[9] = getRegion3x4(ECALRegion3x4_3, 15);
        sort_clusterIn[10] = getRegion3x4(ECALRegion3x4_3, 15);
        sort_clusterIn[11] = getRegion3x4(ECALRegion3x4_3, 15);

// stitch clusters

        stitchClusters(sort_clusterIn, sort_clusterInStitched);

// sorting of clusters, we have now 12 clusters, highest 6
// will be sent to SLRB, they come last from
// sorter

        bitonicSort16(sort_clusterInStitched, sort_clusterOut);
//for(loop i; i<16; i++) cout << sort_clusterOut[i].clusterEnergy() << endl ;
// build ECAL towers with unclustered energy
// keep only 6 highest clusters, others return back to towers

        ap_uint<12> towerEtECAL1[12];
        ap_uint<12> towerEtECAL2[12];
        ap_uint<12> towerEtECAL3[12];
        #pragma HLS ARRAY_PARTITION variable=towerEtECAL1 complete dim=0
        #pragma HLS ARRAY_PARTITION variable=towerEtECAL2 complete dim=0
        #pragma HLS ARRAY_PARTITION variable=towerEtECAL3 complete dim=0

        getTowerEt(ECALRegion3x4_1, towerEtECAL1);
        getTowerEt(ECALRegion3x4_2, towerEtECAL2);
        getTowerEt(ECALRegion3x4_3, towerEtECAL3);


        ap_uint<12> SumE1 ;
        ap_uint<12> SumE2 ;
        ap_uint<12> SumE3 ;

        ap_uint<12> SumE1c[6] ;
        ap_uint<12> SumE2c[6] ;
        ap_uint<12> SumE3c[6] ;
        #pragma HLS ARRAY_PARTITION variable=SumE1c complete dim=0
        #pragma HLS ARRAY_PARTITION variable=SumE2c complete dim=0
        #pragma HLS ARRAY_PARTITION variable=SumE3c complete dim=0

/* -9 takes into account eta offset of the clusters in this SLR */

        for(loop i=0; i<12; i++){
//1                #pragma HLS UNROLL
         for(loop k=4; k<10; k++){
                #pragma HLS UNROLL
          if(sort_clusterOut[k].clusterEnergy() > 0 && (sort_clusterOut[k].towerEta()-9)*4+sort_clusterOut[k].towerPhi() == i){
          SumE1c[k-4] = sort_clusterOut[k].clusterEnergy() ;
          sort_clusterOut[k] = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0) ;
          }
          else SumE1c[k-4] = 0;
          if(sort_clusterOut[k].clusterEnergy() > 0 && (sort_clusterOut[k].towerEta()-9)*4+sort_clusterOut[k].towerPhi() == i+12){
          SumE2c[k-4] = sort_clusterOut[k].clusterEnergy() ;
          sort_clusterOut[k] = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0) ;
          }
          else SumE2c[k-4] = 0;
          if(sort_clusterOut[k].clusterEnergy() > 0 && (sort_clusterOut[k].towerEta()-9)*4+sort_clusterOut[k].towerPhi() == i+24){
          SumE3c[k-4] = sort_clusterOut[k].clusterEnergy() ;
          sort_clusterOut[k] = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0) ;
          }
          else SumE3c[k-4] = 0;
         }
         SumE1 = towerEtECAL1[i] + SumE1c[0] + SumE1c[1]+ SumE1c[2]+ SumE1c[3]+ SumE1c[4]+ SumE1c[5] ;
         SumE2 = towerEtECAL2[i] + SumE2c[0] + SumE2c[1]+ SumE2c[2]+ SumE2c[3]+ SumE2c[4]+ SumE2c[5] ;
         SumE3 = towerEtECAL3[i] + SumE3c[0] + SumE3c[1]+ SumE3c[2]+ SumE3c[3]+ SumE3c[4]+ SumE3c[5] ;
          towerEt[i] = tower_t(SumE1, 0, 0);
          towerEt[i+12] = tower_t(SumE2, 0, 0);
          towerEt[i+24] = tower_t(SumE3, 0, 0);
        }



                        Cluster inClusterLink0[3];
                Cluster inClusterLink1[3];

                tower_t inTowerLink0[16];
                tower_t inTowerLink1[16];

                #pragma HLS ARRAY_PARTITION variable=inClusterLink0 complete dim=0
                #pragma HLS ARRAY_PARTITION variable=inClusterLink0 complete dim=0
                #pragma HLS ARRAY_PARTITION variable=inTowerLink0 complete dim=0
                #pragma HLS ARRAY_PARTITION variable=inTowerLink1 complete dim=0

                for(loop oLink=0; oLink<16; oLink++){
                        #pragma HLS unroll
                        inTowerLink0[oLink] = towerEt[oLink];
                        inTowerLink1[oLink] = towerEt[oLink+16];
                }

                for(loop oLink=0; oLink<3; oLink++){
                        #pragma HLS unroll
                        inClusterLink0[oLink] = sort_clusterOut[oLink+10];
                        inClusterLink1[oLink] = sort_clusterOut[oLink+13];
                }

                                link_out[0] = 0;
                                link_out[1] = 0;

                /*---------------------------------link 1------------------------------------*/
                 processOutLink(inClusterLink0, inTowerLink0, link_out[0]);

                /*---------------------------------link 2------------------------------------*/
                 processOutLink(inClusterLink1, inTowerLink1, link_out[1]);

}
