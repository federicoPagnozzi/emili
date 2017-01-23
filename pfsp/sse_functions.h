//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#include <iostream>
#include <xmmintrin.h>
#include <emmintrin.h>

//#define ENABLE_SSE
/** head and tail computation for taillard's acceleration*/
inline void computeHEADandTAIL(std::vector<int> &sol,std::vector< std::vector < int > >& head,std::vector< std::vector < int > >& tail,const std::vector<std::vector< long> >& processingTimesMatrix,int nbJob, int nbMac)
{
    /** Permutation flowshop Tail and Head matrices computation using SSE instructions
     **/
    /** Each sse register can contain 4 float so the computation is divided in groups of 4 jobs*/
    int r4 = nbJob%4;
    int lambda_number =  nbJob/4; // if ( nbjob%4==0) lambda_number = nbjob/4 else lambda_number = nbjob/4+1;

    int j=1;
    int tj = nbJob;
    std::vector<float> L(nbMac+1,0); // the makespan for each machine of the fourth job in the last group ( at the beginning is zero)
    std::vector<float> tL(nbMac+1,0);
    float res[4] __attribute__((aligned(16)));
    int* k = (int*)res;
        k[0] = 0xffffffff;
        k[1] = 0xffffffff;
        k[2] = 0xffffffff;
        k[3] = 0;
    __m128 mask1110 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0xffffffff;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1100 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1000 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0xffffffff;
    __m128 mask0111 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0110 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0010 = _mm_load_ps(res);
    __m128 K = _mm_setzero_ps();
    __m128 tK = _mm_setzero_ps();
    for(int l = 0 ; l < lambda_number ; l++)
    {
        /** At the beginning
         * J1   J2   J3   J4
        K  T1,1 T2,1 T3,1 T4,1
        L2 T1,2 T2,2 T3,2 T4,2
        L3 T1,3 T2,3 T3,3 T4,3
        .. ..   ..   ..   ..
        LM T1,M T2,M T3,M T4,M

        res[ 0 , 0 , 0 , 0]

        */
        //Initializations
        int j1 = sol[j],j2=sol[j+1],j3=sol[j+2],j4=sol[j+3];
        int tj1 = sol[tj],tj2=sol[tj-1],tj3=sol[tj-2],tj4=sol[tj-3];
        __m128 makespan, mc,mcw;
        __m128 tailspan,tc,tcw;
        makespan = _mm_set_ps(processingTimesMatrix[j4][1],
                processingTimesMatrix[j3][1],
                processingTimesMatrix[j2][1],
                processingTimesMatrix[j1][1]); // load the values in the registers
        mc = makespan; // copy the value in another register

        /** First machine HEAD
         *
         * */
        // first add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a3,a0,a1,a2
        mc = _mm_and_ps(mc,mask0111);         // 0,a0,a1,a2
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1,a2+a3
        // Second add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a2,0,a0,a1
        mc = _mm_and_ps(mc,mask0111);         // 0,0,a0,a1
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1
        // Third add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a1,0,0,a0
        mc = _mm_and_ps(mc,mask0111);         // 0,0,0,a0
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1+a0
        makespan = _mm_add_ps(makespan,K);    // a0+k,a0+a1+k,a2+a1+a0+k,a2+a3+a1+a0+k

        K = _mm_shuffle_ps(makespan,makespan,0xFF);
        _mm_store_ps(res,makespan);

        head[1][j] = res[0];
        head[1][j+1] = res[1];
        head[1][j+2] = res[2];
        head[1][j+3] = res[3];
        /*First machine TAIL
         *
         * */

        tailspan = _mm_set_ps(processingTimesMatrix[tj4][nbMac],
                processingTimesMatrix[tj3][nbMac],
                processingTimesMatrix[tj2][nbMac],
                processingTimesMatrix[tj1][nbMac]);
        // load the values in the registers
        tc = tailspan; // copy the value in another register
                                              // a3,a2,a1,a0
        // first add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a0,a3,a2,a1
        tc = _mm_and_ps(tc,mask0111);         // 0,a3,a2,a1
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1,a0+a1
        // Second add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a1,0,a3,a2
        tc = _mm_and_ps(tc,mask0111);         // 0,0,a3,a2
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1+a3,a2+a0+a1
        // Third add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a2,0,0,a3
        tc = _mm_and_ps(tc,mask0111);         // 0,0,0,a3
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a1+a2+a3,a1+a0+a2+a3
        tailspan = _mm_add_ps(tailspan,tK);    // a3+k,a3+a2+k,a1+a2+a3+k,a1+a0+a2+a3+k

        tK = _mm_shuffle_ps(tailspan,tailspan,0xFF);
        _mm_store_ps(res,tailspan);

        tail[nbMac][tj] = res[0];
        tail[nbMac][tj-1] = res[1];
        tail[nbMac][tj-2] = res[2];
        tail[nbMac][tj-3] = res[3];

        /*The other machines
         *
         **/

        /* Filling the pipeline HEAD
         *
         * */
        int m=2;
                                                                //makespan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        mcw = _mm_set_ps(0,0,0,L[2]);                           // mcw -> [L2 ,0,0,0]
        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,0,processingTimesMatrix[j1][m]);                    // mcw -> [T1,2,0,0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,makespan);
        head[2][j] = res[0];

        // second row
        mcw = _mm_set_ps(0,0,0,L[3]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1000);                           // mc -> [C1,2 , 0   , 0 , 0]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0   , C1,2, 0 , 0]
        mcw = _mm_add_ps(mcw,mc);                               // mcw ->[L3   , C1,2, 0 , 0]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,processingTimesMatrix[j2][m],
                         processingTimesMatrix[j1][m+1]);       // mcw -> [ T1,3, T2,2 , 0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,makespan);
        head[3][j] = res[0];
        head[2][j+1] = res[1];
        // third row
        mcw = _mm_set_ps(0,0,0,L[4]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1100);                           // mc -> [C1,3 , C2,2, 0, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,3, C2,2, 0 ]
        mcw = _mm_add_ps(mcw,mc);                               // mcw -> [L4, C1,3, C2,2  , 0 ]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        mcw = _mm_set_ps(0,processingTimesMatrix[j3][m],
                           processingTimesMatrix[j2][m+1],
                processingTimesMatrix[j1][m+2]);                // mcw -> [ T1,4, T2,3 , T3,2,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,makespan);
        head[4][j] = res[0];
        head[3][j+1] = res[1];
        head[2][j+2] = res[2];
        /* Filling the pipeline TAIL
         *
         * */
        int tm = nbMac-1;
                                                                //tailspan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        tcw = _mm_set_ps(0,0,0,tL[tm]);                           // tcw -> [L2 ,0,0,0]
        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,0,processingTimesMatrix[tj1][tm]);  // tcw -> [T1,2,0,0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm][tj] = res[0];

        // second row
        tcw = _mm_set_ps(0,0,0,tL[tm-1]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1000);                           // tc -> [C1,2 , 0   , 0 , 0]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0   , C1,2, 0 , 0]
        tcw = _mm_add_ps(tcw,tc);                               // tcw ->[L3   , C1,2, 0 , 0]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,processingTimesMatrix[tj2][tm],
                         processingTimesMatrix[tj1][tm-1]);       // tcw -> [ T1,3, T2,2 , 0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,tailspan);
        tail[tm-1][tj] = res[0];
        tail[tm][tj-1] = res[1];
        // third row
        tcw = _mm_set_ps(0,0,0,tL[tm-2]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1100);                           // tc -> [C1,3 , C2,2, 0, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,3, C2,2, 0 ]
        tcw = _mm_add_ps(tcw,tc);                               // tcw -> [L4, C1,3, C2,2  , 0 ]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        tcw = _mm_set_ps(0,processingTimesMatrix[tj3][tm],
                           processingTimesMatrix[tj2][tm-1],
                processingTimesMatrix[tj1][tm-2]);                // tcw -> [ T1,4, T2,3 , T3,2,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj] = res[0];
        tail[tm-1][tj-1] = res[1];
        tail[tm][tj-2] = res[2];
        //other rows
        tm = tm-3;
        for(m = 5; m <= nbMac ; m++)
        {
            /*HEAD
             *
             * */
            // m row
            mcw = _mm_set_ps(0,0,0,L[m]);                       // setup vec for compares
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                       // mc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                    // mc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            mcw = _mm_add_ps(mcw,mc);                           // mcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            makespan = _mm_max_ps(mcw,makespan);                // makespan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            mcw = _mm_set_ps(processingTimesMatrix[j4][m-3],
                    processingTimesMatrix[j3][m-2],
                               processingTimesMatrix[j2][m-1],
                    processingTimesMatrix[j1][m]);              // mcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            makespan = _mm_add_ps(makespan,mcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,makespan);
            head[m][j] = res[0];
            head[m-1][j+1] = res[1];
            head[m-2][j+2] = res[2];
            head[m-3][j+3] = res[3];
            L[m-3] = res[3];

            /*TAIL
             *
             * */

            // m row
            tcw = _mm_set_ps(0,0,0,tL[tm]);                       // setup vec for compares
            tc = tailspan;
            tc = _mm_and_ps(tc,mask1110);                       // tc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            tc = _mm_shuffle_ps(tc,tc,0x93);                    // tc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            tcw = _mm_add_ps(tcw,tc);                           // tcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            tailspan = _mm_max_ps(tcw,tailspan);                // tailspan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            tcw = _mm_set_ps(processingTimesMatrix[tj4][tm+3],
                    processingTimesMatrix[tj3][tm+2],
                               processingTimesMatrix[tj2][tm+1],
                    processingTimesMatrix[tj1][tm]);              // tcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            tailspan = _mm_add_ps(tailspan,tcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,tailspan);
            tail[tm][tj] = res[0];
            tail[tm+1][tj-1] = res[1];
            tail[tm+2][tj-2] = res[2];
            tail[tm+3][tj-3] = res[3];
            tL[tm+3] = res[3];
            tm--;
        }
        /*
         * HEAD finishing...
         * */
                                                        //makespan -> [C1,M , C2,M-1, C3,M-2, C4,M-3]
        // m - 3
        m = nbMac;
        mc = makespan;
        mc = _mm_and_ps(mc,mask1110);                           // mc -> [C1,M , C2,M-1, C3,M-2, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,M , C2,M-1, C3,M-2]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m-2],
                           processingTimesMatrix[j3][m-1],
                processingTimesMatrix[j2][m],0);                // mcw -> [ 0, T2,M , T3,M-1,T4,M-2]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
                                                                //makespan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
        _mm_store_ps(res,makespan);
        head[m][j+1] = res[1];
        head[m-1][j+2] = res[2];
        head[m-2][j+3] = res[3];
        L[m-2] = res[3];
        // m - 2
        mc = makespan;
        mc = _mm_and_ps(mc,mask0110);                           // mc -> [0 , C2,M, C3,M-1, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , C2,M, C3,M-1]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m-1],
                processingTimesMatrix[j3][m],0,0);              // mcw -> [ 0, 0 , T3,M,T4,M-1]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M-1]
        _mm_store_ps(res,makespan);
        head[m][j+2] = res[2];
        head[m-1][j+3] = res[3];
        L[m-1] = res[3];
        //m - 1
        mc = makespan;
        mc = _mm_and_ps(mc,mask0010);                           // mc -> [0 , 0, C3,M, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , 0, C3,M]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , C3,M , max(C3,M, C4,M-1) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m],
                         0,0,0);                                // mcw -> [ 0, 0 , 0,T4,M]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M]
        _mm_store_ps(res,makespan);
        head[m][j+3] = res[3];
        L[m] = res[3];
        j+=4;

        /*TAIL finishing
         *
         * */
        // m - 3
        tm = 3;
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1110);                           // tc -> [C1,M , C2,M-1, C3,M-2, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,M , C2,M-1, C3,M-2]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm],
                           processingTimesMatrix[tj3][tm-1],
                processingTimesMatrix[tj2][tm-2],0);                // tcw -> [ 0, T2,M , T3,M-1,T4,M-2]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
                                                                //tailspan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-1] = res[1];
        tail[tm-1][tj-2] = res[2];
        tail[tm][tj-3] = res[3];
        tL[tm] = res[3];
        // m - 2
        tc = tailspan;
        tc = _mm_and_ps(tc,mask0110);                           // tc -> [0 , C2,M, C3,M-1, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , C2,M, C3,M-1]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm-1],
                processingTimesMatrix[tj3][tm-2],0,0);              // tcw -> [ 0, 0 , T3,M,T4,M-1]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,M, C2,M, C3,M, C4,M-1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-2] = res[2];
        tail[tm-1][tj-3] = res[3];
        tL[tm-1] = res[3];
        //m - 1
        tc = tailspan;
        tc = _mm_and_ps(tc,mask0010);                           // tc -> [0 , 0, C3,M, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , 0, C3,M]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , C3,M , max(C3,M, C4,M-1) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm-2],
                         0,0,0);                                // tcw -> [ 0, 0 , 0,T4,M]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,M, C2,M, C3,M, C4,M]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-3] = res[3];
        tL[tm-2] = res[3];
        tj-=4;
        /* At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }
// What if the number of jobs cannot be divide by 4?
    if(r4 > 0)
    {
        //Initializations
        int j1 = sol[j],j2= r4>=2?sol[j+1]:0,j3=r4>=3?sol[j+2]:0;
        int tj1 = sol[tj],tj2=r4>=2?sol[tj-1]:0,tj3=r4>=3?sol[tj-2]:0;
        __m128 makespan, mc,mcw;
        __m128 tailspan,tc,tcw;
        makespan = _mm_set_ps(0,
                processingTimesMatrix[j3][1],
                processingTimesMatrix[j2][1],
                processingTimesMatrix[j1][1]); // load the values in the registers
        mc = makespan; // copy the value in another register

        /*First machine HEAD
         *
         * */
        // first add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a3,a0,a1,a2
        mc = _mm_and_ps(mc,mask0111);         // 0,a0,a1,a2
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1,a2+a3
        // Second add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a2,0,a0,a1
        mc = _mm_and_ps(mc,mask0111);         // 0,0,a0,a1
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1

        makespan = _mm_add_ps(makespan,K);    // a0+k,a0+a1+k,a2+a1+a0+k,a2+a3+a1+a0+k


        _mm_store_ps(res,makespan);

        head[1][j] = res[0];
        if(r4 > 1)
        {
            head[1][j+1] = res[1];
            if(r4 > 2)
            head[1][j+2] = res[2];
        }
        /*First machine TAIL
         *
         * */

        tailspan = _mm_set_ps(0,
                processingTimesMatrix[tj3][nbMac],
                processingTimesMatrix[tj2][nbMac],
                processingTimesMatrix[tj1][nbMac]); // load the values in the registers
        tc = tailspan; // copy the value in another register
                                              // a3,a2,a1,a0
        // first add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a0,a3,a2,a1
        tc = _mm_and_ps(tc,mask0111);         // 0,a3,a2,a1
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1,a0+a1
        // Second add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a1,0,a3,a2
        tc = _mm_and_ps(tc,mask0111);         // 0,0,a3,a2
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1+a3,a2+a0+a1

        tailspan = _mm_add_ps(tailspan,tK);    // a3+k,a3+a2+k,a1+a2+a3+k,a1+a0+a2+a3+k

        _mm_store_ps(res,tailspan);

        tail[nbMac][tj] = res[0];
        if(r4 > 1)
        {
            tail[nbMac][tj-1] = res[1];
            if(r4 > 2)
            tail[nbMac][tj-2] = res[2];
        }

        /*The other machines
         *
         **/

        /* Filling the pipeline HEAD
         *
         * */
        int m=2;
                                                                //makespan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        mcw = _mm_set_ps(0,0,0,L[2]);                           // mcw -> [L2 ,0,0,0]
        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,0,processingTimesMatrix[j1][m]);                    // mcw -> [T1,2,0,0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,makespan);
        head[2][j] = res[0];

        // second row
        mcw = _mm_set_ps(0,0,0,L[3]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1000);                           // mc -> [C1,2 , 0   , 0 , 0]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0   , C1,2, 0 , 0]
        mcw = _mm_add_ps(mcw,mc);                               // mcw ->[L3   , C1,2, 0 , 0]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,processingTimesMatrix[j2][m],
                         processingTimesMatrix[j1][m+1]);       // mcw -> [ T1,3, T2,2 , 0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,makespan);
        head[3][j] = res[0];
        if(r4>1)
        head[2][j+1] = res[1];
        // third row
        mcw = _mm_set_ps(0,0,0,L[4]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1100);                           // mc -> [C1,3 , C2,2,    0, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [   0 , C1,3, C2,2, 0 ]
        mcw = _mm_add_ps(mcw,mc);                               // mcw ->[  L4 , C1,3, C2,2, 0 ]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        mcw = _mm_set_ps(0,processingTimesMatrix[j3][m],
                           processingTimesMatrix[j2][m+1],
                processingTimesMatrix[j1][m+2]);                // mcw -> [ T1,4, T2,3 , T3,2,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,makespan);
        head[4][j] = res[0];
        if(r4>1)
        head[3][j+1] = res[1];
        if(r4>2)
        head[2][j+2] = res[2];
        /* Filling the pipeline TAIL
         *
         * */
        int tm = nbMac-1;
                                                                //tailspan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        tcw = _mm_set_ps(0,0,0,tL[tm]);                           // tcw -> [L2 ,0,0,0]
        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,0,processingTimesMatrix[tj1][tm]);  // tcw -> [T1,2,0,0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm][tj] = res[0];

        // second row
        tcw = _mm_set_ps(0,0,0,tL[tm-1]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1000);                           // tc -> [C1,2 , 0   , 0 , 0]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0   , C1,2, 0 , 0]
        tcw = _mm_add_ps(tcw,tc);                               // tcw ->[L3   , C1,2, 0 , 0]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,processingTimesMatrix[tj2][tm],
                         processingTimesMatrix[tj1][tm-1]);       // tcw -> [ T1,3, T2,2 , 0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,tailspan);
        tail[tm-1][tj] = res[0];
        if(r4>1)
        tail[tm][tj-1] = res[1];
        // third row
        tcw = _mm_set_ps(0,0,0,tL[tm-2]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1100);                           // tc -> [C1,3 , C2,2, 0, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,3, C2,2, 0 ]
        tcw = _mm_add_ps(tcw,tc);                               // tcw -> [L4, C1,3, C2,2  , 0 ]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        tcw = _mm_set_ps(0,processingTimesMatrix[tj3][tm],
                           processingTimesMatrix[tj2][tm-1],
                processingTimesMatrix[tj1][tm-2]);                // tcw -> [ T1,4, T2,3 , T3,2,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj] = res[0];
        if(r4>1)
        tail[tm-1][tj-1] = res[1];
        if(r4>2)
        tail[tm][tj-2] = res[2];
        //other rows
        tm = tm-3;
        for(m = 5; m <= nbMac ; m++)
        {
            /*HEAD
             *
             * */
            // m row
            mcw = _mm_set_ps(0,0,0,L[m]);                       // setup vec for compares
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                       // mc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                    // mc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            mcw = _mm_add_ps(mcw,mc);                           // mcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            makespan = _mm_max_ps(mcw,makespan);                // makespan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            mcw = _mm_set_ps(0,
                    processingTimesMatrix[j3][m-2],
                               processingTimesMatrix[j2][m-1],
                    processingTimesMatrix[j1][m]);              // mcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            makespan = _mm_add_ps(makespan,mcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,makespan);
            head[m][j] = res[0];
            if(r4>1)
            head[m-1][j+1] = res[1];
            if(r4>2)
            head[m-2][j+2] = res[2];
            /*TAIL
             *
             * */
            tcw = _mm_set_ps(0,0,0,tL[tm]);                       // setup vec for compares
            tc = tailspan;
            tc = _mm_and_ps(tc,mask1110);                       // tc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            tc = _mm_shuffle_ps(tc,tc,0x93);                    // tc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            tcw = _mm_add_ps(tcw,tc);                           // tcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            tailspan = _mm_max_ps(tcw,tailspan);                // tailspan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            tcw = _mm_set_ps(0,
                    processingTimesMatrix[tj3][tm+2],
                               processingTimesMatrix[tj2][tm+1],
                    processingTimesMatrix[tj1][tm]);              // tcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            tailspan = _mm_add_ps(tailspan,tcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,tailspan);
            tail[tm][tj] = res[0];
            if(r4>1)
            tail[tm+1][tj-1] = res[1];
            if(r4>2)
            tail[tm+2][tj-2] = res[2];
            tm--;
        }
        /*
         * HEAD finishing...
         * */
                                                        //makespan -> [C1,M , C2,M-1, C3,M-2, C4,M-3]
        // m - 3
        if(r4>1)
        {
            m = nbMac;
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                           // mc -> [C1,M , C2,M-1, C3,M-2, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,M , C2,M-1, C3,M-2]
            makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
            mcw = _mm_set_ps(0,
                             processingTimesMatrix[j3][m-1],
                    processingTimesMatrix[j2][m],0);                // mcw -> [ 0, T2,M , T3,M-1,T4,M-2]
            makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
            //makespan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
            _mm_store_ps(res,makespan);
            head[m][j+1] = res[1];
            if(r4>2)
            {
                head[m-1][j+2] = res[2];
                // m - 2
                mc = makespan;
                mc = _mm_and_ps(mc,mask0110);                           // mc -> [0 , C2,M, C3,M-1, 0 ]
                mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , C2,M, C3,M-1]
                makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
                mcw = _mm_set_ps(0,
                        processingTimesMatrix[j3][m],0,0);              // mcw -> [ 0, 0 , T3,M,T4,M-1]
                makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                //makespan -> [ C1,M, C2,M, C3,M, C4,M-1]
                _mm_store_ps(res,makespan);
                head[m][j+2] = res[2];
            }
        }
        /*TAIL finishing
         *
         * */
        // m - 3
        if(r4 > 1)
        {
            tm = 3;
            tc = tailspan;
            tc = _mm_and_ps(tc,mask1110);                           // tc -> [C1,M , C2,M-1, C3,M-2, 0 ]
            tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,M , C2,M-1, C3,M-2]
            tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
            tcw = _mm_set_ps(0,
                             processingTimesMatrix[tj3][tm-1],
                    processingTimesMatrix[tj2][tm-2],0);                // tcw -> [ 0, T2,M , T3,M-1,T4,M-2]
            tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
            //tailspan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
            _mm_store_ps(res,tailspan);
            tail[tm-2][tj-1] = res[1];
            if(r4>2)
            {
                tail[tm-1][tj-2] = res[2];
                // m - 2
                tc = tailspan;
                tc = _mm_and_ps(tc,mask0110);                           // tc -> [0 , C2,M, C3,M-1, 0 ]
                tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , C2,M, C3,M-1]
                tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
                tcw = _mm_set_ps(0,
                        processingTimesMatrix[tj3][tm-2],0,0);              // tcw -> [ 0, 0 , T3,M,T4,M-1]
                tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                //tailspan -> [ C1,M, C2,M, C3,M, C4,M-1]
                _mm_store_ps(res,tailspan);
                tail[tm-2][tj-2] = res[2];
            }
        }
        /* At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }

}
// this function calculates only the head
inline void computeHEAD(std::vector<int> &sol,std::vector< std::vector < int > >& head,const std::vector<std::vector< long> >& processingTimesMatrix,int nbJob, int nbMac)
{
    /* Permutation flowshop Head matrix computation using SSE instructions
     **/
    // Each sse register can contain 4 float so the computation is divided in groups of 4 jobs
    int r4 = nbJob%4;
    int lambda_number =  nbJob/4; // if ( nbjob%4==0) lambda_number = nbjob/4 else lambda_number = nbjob/4+1;

    int j=1;
    std::vector<float> L(nbMac+1,0); // the makespan for each machine of the fourth job in the last group ( at the beginning is zero)
    float res[4] __attribute__((aligned(16)));
    int* k = (int*)res;
        k[0] = 0xffffffff;
        k[1] = 0xffffffff;
        k[2] = 0xffffffff;
        k[3] = 0;
    __m128 mask1110 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0xffffffff;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1100 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1000 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0xffffffff;
    __m128 mask0111 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0110 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0010 = _mm_load_ps(res);
    __m128 K = _mm_setzero_ps();

    for(int l = 0 ; l < lambda_number ; l++)
    {
        /* At the beginning
         * J1   J2   J3   J4
        K  T1,1 T2,1 T3,1 T4,1
        L2 T1,2 T2,2 T3,2 T4,2
        L3 T1,3 T2,3 T3,3 T4,3
        .. ..   ..   ..   ..
        LM T1,M T2,M T3,M T4,M

        res[ 0 , 0 , 0 , 0]

        */
        //Initializations
        int j1 = sol[j],j2=sol[j+1],j3=sol[j+2],j4=sol[j+3];

        __m128 makespan, mc,mcw;

        makespan = _mm_set_ps(processingTimesMatrix[j4][1],
                processingTimesMatrix[j3][1],
                processingTimesMatrix[j2][1],
                processingTimesMatrix[j1][1]); // load the values in the registers
        mc = makespan; // copy the value in another register

        /*First machine HEAD
         *
         * */
        // first add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a3,a0,a1,a2
        mc = _mm_and_ps(mc,mask0111);         // 0,a0,a1,a2
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1,a2+a3
        // Second add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a2,0,a0,a1
        mc = _mm_and_ps(mc,mask0111);         // 0,0,a0,a1
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1
        // Third add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a1,0,0,a0
        mc = _mm_and_ps(mc,mask0111);         // 0,0,0,a0
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1+a0
        makespan = _mm_add_ps(makespan,K);    // a0+k,a0+a1+k,a2+a1+a0+k,a2+a3+a1+a0+k

        K = _mm_shuffle_ps(makespan,makespan,0xFF);
        _mm_store_ps(res,makespan);

        head[1][j] = res[0];
        head[1][j+1] = res[1];
        head[1][j+2] = res[2];
        head[1][j+3] = res[3];

        /*The other machines
         *
         **/

        /* Filling the pipeline HEAD
         *
         * */
        int m=2;
                                                                //makespan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        mcw = _mm_set_ps(0,0,0,L[2]);                           // mcw -> [L2 ,0,0,0]
        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,0,processingTimesMatrix[j1][m]);                    // mcw -> [T1,2,0,0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,makespan);
        head[2][j] = res[0];

        // second row
        mcw = _mm_set_ps(0,0,0,L[3]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1000);                           // mc -> [C1,2 , 0   , 0 , 0]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0   , C1,2, 0 , 0]
        mcw = _mm_add_ps(mcw,mc);                               // mcw ->[L3   , C1,2, 0 , 0]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,processingTimesMatrix[j2][m],
                         processingTimesMatrix[j1][m+1]);       // mcw -> [ T1,3, T2,2 , 0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,makespan);
        head[3][j] = res[0];
        head[2][j+1] = res[1];
        // third row
        mcw = _mm_set_ps(0,0,0,L[4]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1100);                           // mc -> [C1,3 , C2,2, 0, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,3, C2,2, 0 ]
        mcw = _mm_add_ps(mcw,mc);                               // mcw -> [L4, C1,3, C2,2  , 0 ]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        mcw = _mm_set_ps(0,processingTimesMatrix[j3][m],
                           processingTimesMatrix[j2][m+1],
                processingTimesMatrix[j1][m+2]);                // mcw -> [ T1,4, T2,3 , T3,2,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,makespan);
        head[4][j] = res[0];
        head[3][j+1] = res[1];
        head[2][j+2] = res[2];

        for(m = 5; m <= nbMac ; m++)
        {
            /*HEAD
             *
             * */
            // m row
            mcw = _mm_set_ps(0,0,0,L[m]);                       // setup vec for compares
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                       // mc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                    // mc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            mcw = _mm_add_ps(mcw,mc);                           // mcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            makespan = _mm_max_ps(mcw,makespan);                // makespan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            mcw = _mm_set_ps(processingTimesMatrix[j4][m-3],
                    processingTimesMatrix[j3][m-2],
                               processingTimesMatrix[j2][m-1],
                    processingTimesMatrix[j1][m]);              // mcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            makespan = _mm_add_ps(makespan,mcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,makespan);
            head[m][j] = res[0];
            head[m-1][j+1] = res[1];
            head[m-2][j+2] = res[2];
            head[m-3][j+3] = res[3];
            L[m-3] = res[3];
        }
        /*
         * HEAD finishing...
         * */
                                                        //makespan -> [C1,M , C2,M-1, C3,M-2, C4,M-3]
        // m - 3
        m = nbMac;
        mc = makespan;
        mc = _mm_and_ps(mc,mask1110);                           // mc -> [C1,M , C2,M-1, C3,M-2, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,M , C2,M-1, C3,M-2]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m-2],
                           processingTimesMatrix[j3][m-1],
                processingTimesMatrix[j2][m],0);                // mcw -> [ 0, T2,M , T3,M-1,T4,M-2]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
                                                                //makespan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
        _mm_store_ps(res,makespan);
        head[m][j+1] = res[1];
        head[m-1][j+2] = res[2];
        head[m-2][j+3] = res[3];
        L[m-2] = res[3];
        // m - 2
        mc = makespan;
        mc = _mm_and_ps(mc,mask0110);                           // mc -> [0 , C2,M, C3,M-1, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , C2,M, C3,M-1]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m-1],
                processingTimesMatrix[j3][m],0,0);              // mcw -> [ 0, 0 , T3,M,T4,M-1]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M-1]
        _mm_store_ps(res,makespan);
        head[m][j+2] = res[2];
        head[m-1][j+3] = res[3];
        L[m-1] = res[3];
        //m - 1
        mc = makespan;
        mc = _mm_and_ps(mc,mask0010);                           // mc -> [0 , 0, C3,M, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , 0, C3,M]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , C3,M , max(C3,M, C4,M-1) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m],
                         0,0,0);                                // mcw -> [ 0, 0 , 0,T4,M]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M]
        _mm_store_ps(res,makespan);
        head[m][j+3] = res[3];
        L[m] = res[3];
        j+=4;

        /* At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }
// What if the number of jobs cannot be divide by 4?
    if(r4 > 0)
    {
        //Initializations
        int j1 = sol[j],j2= r4>=2?sol[j+1]:0,j3=r4>=3?sol[j+2]:0;
        __m128 makespan, mc,mcw;
        makespan = _mm_set_ps(0,
                processingTimesMatrix[j3][1],
                processingTimesMatrix[j2][1],
                processingTimesMatrix[j1][1]); // load the values in the registers
        mc = makespan; // copy the value in another register

        /*First machine HEAD
         *
         * */
        // first add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a3,a0,a1,a2
        mc = _mm_and_ps(mc,mask0111);         // 0,a0,a1,a2
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1,a2+a3
        // Second add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a2,0,a0,a1
        mc = _mm_and_ps(mc,mask0111);         // 0,0,a0,a1
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1

        makespan = _mm_add_ps(makespan,K);    // a0+k,a0+a1+k,a2+a1+a0+k,a2+a3+a1+a0+k


        _mm_store_ps(res,makespan);

        head[1][j] = res[0];
        if(r4 > 1)
        {
            head[1][j+1] = res[1];
            if(r4 > 2)
            head[1][j+2] = res[2];
        }

        /*The other machines
         *
         **/

        /* Filling the pipeline HEAD
         *
         * */
        int m=2;
                                                                //makespan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        mcw = _mm_set_ps(0,0,0,L[2]);                           // mcw -> [L2 ,0,0,0]
        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,0,processingTimesMatrix[j1][m]);                    // mcw -> [T1,2,0,0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,makespan);
        head[2][j] = res[0];

        // second row
        mcw = _mm_set_ps(0,0,0,L[3]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1000);                           // mc -> [C1,2 , 0   , 0 , 0]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0   , C1,2, 0 , 0]
        mcw = _mm_add_ps(mcw,mc);                               // mcw ->[L3   , C1,2, 0 , 0]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,processingTimesMatrix[j2][m],
                         processingTimesMatrix[j1][m+1]);       // mcw -> [ T1,3, T2,2 , 0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,makespan);
        head[3][j] = res[0];
        if(r4>1)
        head[2][j+1] = res[1];
        // third row
        mcw = _mm_set_ps(0,0,0,L[4]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1100);                           // mc -> [C1,3 , C2,2, 0, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,3, C2,2, 0 ]
        mcw = _mm_add_ps(mcw,mc);                               // mcw -> [L4, C1,3, C2,2  , 0 ]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        mcw = _mm_set_ps(0,processingTimesMatrix[j3][m],
                           processingTimesMatrix[j2][m+1],
                processingTimesMatrix[j1][m+2]);                // mcw -> [ T1,4, T2,3 , T3,2,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,makespan);
        head[4][j] = res[0];
        if(r4>1)
        head[3][j+1] = res[1];
        if(r4>2)
        head[2][j+2] = res[2];

        for(m = 5; m <= nbMac ; m++)
        {
            /*HEAD
             *
             * */
            // m row
            mcw = _mm_set_ps(0,0,0,L[m]);                       // setup vec for compares
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                       // mc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                    // mc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            mcw = _mm_add_ps(mcw,mc);                           // mcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            makespan = _mm_max_ps(mcw,makespan);                // makespan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            mcw = _mm_set_ps(0,
                    processingTimesMatrix[j3][m-2],
                               processingTimesMatrix[j2][m-1],
                    processingTimesMatrix[j1][m]);              // mcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            makespan = _mm_add_ps(makespan,mcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,makespan);
            head[m][j] = res[0];
            if(r4>1)
            head[m-1][j+1] = res[1];
            if(r4>2)
            head[m-2][j+2] = res[2];
        }
        /*
         * HEAD finishing...
         * */
                                                        //makespan -> [C1,M , C2,M-1, C3,M-2, C4,M-3]
        // m - 3
        if(r4>1)
        {
            m = nbMac;
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                           // mc -> [C1,M , C2,M-1, C3,M-2, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,M , C2,M-1, C3,M-2]
            makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
            mcw = _mm_set_ps(0,
                             processingTimesMatrix[j3][m-1],
                    processingTimesMatrix[j2][m],0);                // mcw -> [ 0, T2,M , T3,M-1,T4,M-2]
            makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
            //makespan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
            _mm_store_ps(res,makespan);
            head[m][j+1] = res[1];
            if(r4>2)
            {
                head[m-1][j+2] = res[2];
                // m - 2
                mc = makespan;
                mc = _mm_and_ps(mc,mask0110);                           // mc -> [0 , C2,M, C3,M-1, 0 ]
                mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , C2,M, C3,M-1]
                makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
                mcw = _mm_set_ps(0,
                        processingTimesMatrix[j3][m],0,0);              // mcw -> [ 0, 0 , T3,M,T4,M-1]
                makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                //makespan -> [ C1,M, C2,M, C3,M, C4,M-1]
                _mm_store_ps(res,makespan);
                head[m][j+2] = res[2];
            }
        }
        /* At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }
}
// this function calculates only the tail
inline void computeTAIL(std::vector<int> &sol,std::vector< std::vector < int > >& tail,const std::vector<std::vector< long> >& processingTimesMatrix,int nbJob, int nbMac)
{
    /* Permutation flowshop Tail matrix computation using SSE instructions
     **/
    // Each sse register can contain 4 float so the computation is divided in groups of 4 jobs
    int r4 = nbJob%4;
    int lambda_number =  nbJob/4; // if ( nbjob%4==0) lambda_number = nbjob/4 else lambda_number = nbjob/4+1;

    int tj = nbJob;
// the makespan for each machine of the fourth job in the last group ( at the beginning is zero)
    std::vector<float> tL(nbMac+1,0);
    float res[4] __attribute__((aligned(16)));
    int* k = (int*)res;
        k[0] = 0xffffffff;
        k[1] = 0xffffffff;
        k[2] = 0xffffffff;
        k[3] = 0;
    __m128 mask1110 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0xffffffff;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1100 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1000 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0xffffffff;
    __m128 mask0111 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0110 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0010 = _mm_load_ps(res);
    __m128 tK = _mm_setzero_ps();
    for(int l = 0 ; l < lambda_number ; l++)
    {
        /* At the beginning
         * J1   J2   J3   J4
        K  T1,1 T2,1 T3,1 T4,1
        L2 T1,2 T2,2 T3,2 T4,2
        L3 T1,3 T2,3 T3,3 T4,3
        .. ..   ..   ..   ..
        LM T1,M T2,M T3,M T4,M

        res[ 0 , 0 , 0 , 0]

        */
        //Initializations
        int tj1 = sol[tj],tj2=sol[tj-1],tj3=sol[tj-2],tj4=sol[tj-3];
        __m128 tailspan,tc,tcw;
        /*First machine TAIL
         *
         * */

        tailspan = _mm_set_ps(processingTimesMatrix[tj4][nbMac],
                processingTimesMatrix[tj3][nbMac],
                processingTimesMatrix[tj2][nbMac],
                processingTimesMatrix[tj1][nbMac]);
        // load the values in the registers
        tc = tailspan; // copy the value in another register
                                              // a3,a2,a1,a0
        // first add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a0,a3,a2,a1
        tc = _mm_and_ps(tc,mask0111);         // 0,a3,a2,a1
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1,a0+a1
        // Second add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a1,0,a3,a2
        tc = _mm_and_ps(tc,mask0111);         // 0,0,a3,a2
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1+a3,a2+a0+a1
        // Third add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a2,0,0,a3
        tc = _mm_and_ps(tc,mask0111);         // 0,0,0,a3
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a1+a2+a3,a1+a0+a2+a3
        tailspan = _mm_add_ps(tailspan,tK);    // a3+k,a3+a2+k,a1+a2+a3+k,a1+a0+a2+a3+k

        tK = _mm_shuffle_ps(tailspan,tailspan,0xFF);
        _mm_store_ps(res,tailspan);

        tail[nbMac][tj] = res[0];
        tail[nbMac][tj-1] = res[1];
        tail[nbMac][tj-2] = res[2];
        tail[nbMac][tj-3] = res[3];

        /*The other machines
         *
         **/
        /* Filling the pipeline TAIL
         *
         * */
        int tm = nbMac-1;
                                                                //tailspan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        tcw = _mm_set_ps(0,0,0,tL[tm]);                           // tcw -> [L2 ,0,0,0]
        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,0,processingTimesMatrix[tj1][tm]);  // tcw -> [T1,2,0,0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm][tj] = res[0];

        // second row
        tcw = _mm_set_ps(0,0,0,tL[tm-1]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1000);                           // tc -> [C1,2 , 0   , 0 , 0]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0   , C1,2, 0 , 0]
        tcw = _mm_add_ps(tcw,tc);                               // tcw ->[L3   , C1,2, 0 , 0]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,processingTimesMatrix[tj2][tm],
                         processingTimesMatrix[tj1][tm-1]);       // tcw -> [ T1,3, T2,2 , 0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,tailspan);
        tail[tm-1][tj] = res[0];
        tail[tm][tj-1] = res[1];
        // third row
        tcw = _mm_set_ps(0,0,0,tL[tm-2]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1100);                           // tc -> [C1,3 , C2,2, 0, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,3, C2,2, 0 ]
        tcw = _mm_add_ps(tcw,tc);                               // tcw -> [L4, C1,3, C2,2  , 0 ]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        tcw = _mm_set_ps(0,processingTimesMatrix[tj3][tm],
                           processingTimesMatrix[tj2][tm-1],
                processingTimesMatrix[tj1][tm-2]);                // tcw -> [ T1,4, T2,3 , T3,2,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj] = res[0];
        tail[tm-1][tj-1] = res[1];
        tail[tm][tj-2] = res[2];
        //other rows

        for(tm = nbMac-4; tm > 0  ; tm--)
        {
            /*TAIL
             *
             * */

            // m row
            tcw = _mm_set_ps(0,0,0,tL[tm]);                       // setup vec for compares
            tc = tailspan;
            tc = _mm_and_ps(tc,mask1110);                       // tc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            tc = _mm_shuffle_ps(tc,tc,0x93);                    // tc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            tcw = _mm_add_ps(tcw,tc);                           // tcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            tailspan = _mm_max_ps(tcw,tailspan);                // tailspan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            tcw = _mm_set_ps(processingTimesMatrix[tj4][tm+3],
                    processingTimesMatrix[tj3][tm+2],
                               processingTimesMatrix[tj2][tm+1],
                    processingTimesMatrix[tj1][tm]);              // tcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            tailspan = _mm_add_ps(tailspan,tcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,tailspan);
            tail[tm][tj] = res[0];
            tail[tm+1][tj-1] = res[1];
            tail[tm+2][tj-2] = res[2];
            tail[tm+3][tj-3] = res[3];
            tL[tm+3] = res[3];

        }

        /*TAIL finishing
         *
         * */
        // m - 3
        tm = 3;
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1110);                           // tc -> [C1,M , C2,M-1, C3,M-2, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,M , C2,M-1, C3,M-2]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm],
                           processingTimesMatrix[tj3][tm-1],
                processingTimesMatrix[tj2][tm-2],0);                // tcw -> [ 0, T2,M , T3,M-1,T4,M-2]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
                                                                //tailspan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-1] = res[1];
        tail[tm-1][tj-2] = res[2];
        tail[tm][tj-3] = res[3];
        tL[tm] = res[3];
        // m - 2
        tc = tailspan;
        tc = _mm_and_ps(tc,mask0110);                           // tc -> [0 , C2,M, C3,M-1, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , C2,M, C3,M-1]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm-1],
                processingTimesMatrix[tj3][tm-2],0,0);              // tcw -> [ 0, 0 , T3,M,T4,M-1]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,M, C2,M, C3,M, C4,M-1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-2] = res[2];
        tail[tm-1][tj-3] = res[3];
        tL[tm-1] = res[3];
        //m - 1
        tc = tailspan;
        tc = _mm_and_ps(tc,mask0010);                           // tc -> [0 , 0, C3,M, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , 0, C3,M]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , C3,M , max(C3,M, C4,M-1) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm-2],
                         0,0,0);                                // tcw -> [ 0, 0 , 0,T4,M]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,M, C2,M, C3,M, C4,M]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-3] = res[3];
        tL[tm-2] = res[3];
        tj-=4;
        /* At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }
// What if the number of jobs cannot be divided by 4?
    if(r4 > 0)
    {
        //Initializations
        int tj1 = sol[tj],tj2=r4>=2?sol[tj-1]:0,tj3=r4>=3?sol[tj-2]:0;
        __m128 tailspan,tc,tcw;
        /*First machine TAIL
         *
         * */

        tailspan = _mm_set_ps(0,
                processingTimesMatrix[tj3][nbMac],
                processingTimesMatrix[tj2][nbMac],
                processingTimesMatrix[tj1][nbMac]); // load the values in the registers
        tc = tailspan; // copy the value in another register
                                              // a3,a2,a1,a0
        // first add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a0,a3,a2,a1
        tc = _mm_and_ps(tc,mask0111);         // 0,a3,a2,a1
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1,a0+a1
        // Second add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a1,0,a3,a2
        tc = _mm_and_ps(tc,mask0111);         // 0,0,a3,a2
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1+a3,a2+a0+a1

        tailspan = _mm_add_ps(tailspan,tK);    // a3+k,a3+a2+k,a1+a2+a3+k,a1+a0+a2+a3+k

        _mm_store_ps(res,tailspan);

        tail[nbMac][tj] = res[0];
        if(r4 > 1)
        {
            tail[nbMac][tj-1] = res[1];
            if(r4 > 2)
            tail[nbMac][tj-2] = res[2];
        }

        /*The other machines
         *
         **/

        /* Filling the pipeline TAIL
         *
         * */
        int tm = nbMac-1;
                                                                //tailspan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        tcw = _mm_set_ps(0,0,0,tL[tm]);                           // tcw -> [L2 ,0,0,0]
        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,0,processingTimesMatrix[tj1][tm]);  // tcw -> [T1,2,0,0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm][tj] = res[0];

        // second row
        tcw = _mm_set_ps(0,0,0,tL[tm-1]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1000);                           // tc -> [C1,2 , 0   , 0 , 0]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0   , C1,2, 0 , 0]
        tcw = _mm_add_ps(tcw,tc);                               // tcw ->[L3   , C1,2, 0 , 0]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,processingTimesMatrix[tj2][tm],
                         processingTimesMatrix[tj1][tm-1]);       // tcw -> [ T1,3, T2,2 , 0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,tailspan);
        tail[tm-1][tj] = res[0];
        if(r4>1)
        tail[tm][tj-1] = res[1];
        // third row
        tcw = _mm_set_ps(0,0,0,tL[tm-2]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1100);                           // tc -> [C1,3 , C2,2, 0, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,3, C2,2, 0 ]
        tcw = _mm_add_ps(tcw,tc);                               // tcw -> [L4, C1,3, C2,2  , 0 ]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        tcw = _mm_set_ps(0,processingTimesMatrix[tj3][tm],
                           processingTimesMatrix[tj2][tm-1],
                processingTimesMatrix[tj1][tm-2]);                // tcw -> [ T1,4, T2,3 , T3,2,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj] = res[0];
        if(r4>1)
        tail[tm-1][tj-1] = res[1];
        if(r4>2)
        tail[tm][tj-2] = res[2];
        //other rows

        for(tm = nbMac-4; tm > 0 ; tm--)
        {
            /*TAIL
             *
             * */
            tcw = _mm_set_ps(0,0,0,tL[tm]);                       // setup vec for compares
            tc = tailspan;
            tc = _mm_and_ps(tc,mask1110);                       // tc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            tc = _mm_shuffle_ps(tc,tc,0x93);                    // tc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            tcw = _mm_add_ps(tcw,tc);                           // tcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            tailspan = _mm_max_ps(tcw,tailspan);                // tailspan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            tcw = _mm_set_ps(0,
                    processingTimesMatrix[tj3][tm+2],
                               processingTimesMatrix[tj2][tm+1],
                    processingTimesMatrix[tj1][tm]);              // tcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            tailspan = _mm_add_ps(tailspan,tcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,tailspan);
            tail[tm][tj] = res[0];
            if(r4>1)
            tail[tm+1][tj-1] = res[1];
            if(r4>2)
            tail[tm+2][tj-2] = res[2];
        }
        /*TAIL finishing
         *
         * */
        // m - 3
        if(r4 > 1)
        {
            tm = 3;
            tc = tailspan;
            tc = _mm_and_ps(tc,mask1110);                           // tc -> [C1,M , C2,M-1, C3,M-2, 0 ]
            tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,M , C2,M-1, C3,M-2]
            tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
            tcw = _mm_set_ps(0,
                             processingTimesMatrix[tj3][tm-1],
                    processingTimesMatrix[tj2][tm-2],0);                // tcw -> [ 0, T2,M , T3,M-1,T4,M-2]
            tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
            //tailspan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
            _mm_store_ps(res,tailspan);
            tail[tm-2][tj-1] = res[1];
            if(r4>2)
            {
                tail[tm-1][tj-2] = res[2];
                // m - 2
                tc = tailspan;
                tc = _mm_and_ps(tc,mask0110);                           // tc -> [0 , C2,M, C3,M-1, 0 ]
                tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , C2,M, C3,M-1]
                tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
                tcw = _mm_set_ps(0,
                        processingTimesMatrix[tj3][tm-2],0,0);              // tcw -> [ 0, 0 , T3,M,T4,M-1]
                tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                //tailspan -> [ C1,M, C2,M, C3,M, C4,M-1]
                _mm_store_ps(res,tailspan);
                tail[tm-2][tj-2] = res[2];
            }
        }
        /* At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }
}
// This function calculates completion times of all jobs
inline void computePMakespans(std::vector<int>& sol,std::vector< long >& previousMachineEndTime,const std::vector<std::vector< long> >& pmat,int nJob, int nbMac,int starting_job,float* L)
{

    /* Permutation flowshop makespan computation using SSE instructions
     **/
    // Each sse register can contain 4 float so the computation is divided in groups of 4 jobs
    int nbJob = nJob-starting_job;
    int r4 = nbJob%4;
    int lambda_number =  r4==0?nbJob/4:(nbJob/4+1); // if ( nbjob%4==0) lambda_number = nbjob/4 else lambda_number = nbjob/4+1;
    if(r4>0)
    {
        for(int i=0;i<(4-r4);i++)
        {
            sol.push_back(0);
            previousMachineEndTime.push_back(0);
        }
    }
    int j=starting_job;
     // the makespan for each machine of the fourth job in the last group ( at the beginning is zero)
    float res[4] __attribute__((aligned(16)));
    int* k = (int*)res;
        k[0] = 0xffffffff;
        k[1] = 0xffffffff;
        k[2] = 0xffffffff;
        k[3] = 0;
    __m128 mask1110 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0xffffffff;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1100 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1000 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0xffffffff;
    __m128 mask0111 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0110 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0010 = _mm_load_ps(res);
    __m128 K = _mm_set_ps1(L[1]);

    for(int l = 0 ; l < lambda_number ; l++)
    {
        /* At the beginning
         * J1   J2   J3   J4
        K  T1,1 T2,1 T3,1 T4,1
        L2 T1,2 T2,2 T3,2 T4,2
        L3 T1,3 T2,3 T3,3 T4,3
        .. ..   ..   ..   ..
        LM T1,M T2,M T3,M T4,M

        res[ 0 , 0 , 0 , 0]

        */
        //Initializations
        int j1 = sol[j],j2=sol[j+1],j3=sol[j+2],j4=sol[j+3];
        __m128 makespan, mc,mcw;
        makespan = _mm_set_ps(pmat[j4][1],pmat[j3][1],pmat[j2][1],pmat[j1][1]); // load the values in the registers
        mc = makespan; // copy the value in another register

        /*First machine
         *
         * */
        // first add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a3,a0,a1,a2
        mc = _mm_and_ps(mc,mask0111);         // 0,a0,a1,a2
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1,a2+a3
        // Second add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a2,0,a0,a1
        mc = _mm_and_ps(mc,mask0111);         // 0,0,a0,a1
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1
        // Third add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a1,0,0,a0
        mc = _mm_and_ps(mc,mask0111);         // 0,0,0,a0
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1+a0
        makespan = _mm_add_ps(makespan,K);    // a0+k,a0+a1+k,a2+a1+a0+k,a2+a3+a1+a0+k

        K = _mm_shuffle_ps(makespan,makespan,0xFF);

        /*The other machines
         *
         **/
        int m=2;
                                                                //makespan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        mcw = _mm_set_ps(0,0,0,L[2]);                           // mcw -> [L2 ,0,0,0]
        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,0,pmat[j1][m]);                    // mcw -> [T1,2,0,0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,2, C2,1 , C3,1 , C4,1]

        // second row
        mcw = _mm_set_ps(0,0,0,L[3]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1000);                           // mc -> [C1,2 , 0   , 0 , 0]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0   , C1,2, 0 , 0]
        mcw = _mm_add_ps(mcw,mc);                               // mcw ->[L3   , C1,2, 0 , 0]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,pmat[j2][m],pmat[j1][m+1]);        // mcw -> [ T1,3, T2,2 , 0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)

        // third row
        mcw = _mm_set_ps(0,0,0,L[4]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1100);                           // mc -> [C1,3 , C2,2, 0, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,3, C2,2, 0 ]
        mcw = _mm_add_ps(mcw,mc);                               // mcw -> [L4, C1,3, C2,2  , 0 ]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        mcw = _mm_set_ps(0,pmat[j3][m],
                           pmat[j2][m+1],pmat[j1][m+2]);        // mcw -> [ T1,4, T2,3 , T3,2,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,4, C2,3, C3,2, C4,1]
        //other rows
        for(m = 5; m <= nbMac ; m++)
        {
            // m row
            mcw = _mm_set_ps(0,0,0,L[m]);                       // setup vec for compares
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                       // mc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                    // mc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            mcw = _mm_add_ps(mcw,mc);                           // mcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            makespan = _mm_max_ps(mcw,makespan);                // makespan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            mcw = _mm_set_ps(pmat[j4][m-3],pmat[j3][m-2],
                               pmat[j2][m-1],pmat[j1][m]);      // mcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            makespan = _mm_add_ps(makespan,mcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,makespan);
            L[m-3] = res[3];
        }
                                                                //makespan -> [C1,M , C2,M-1, C3,M-2, C4,M-3]
        // m - 3
        m = nbMac;
        mc = makespan;
        mc = _mm_and_ps(mc,mask1110);                           // mc -> [C1,M , C2,M-1, C3,M-2, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,M , C2,M-1, C3,M-2]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
        mcw = _mm_set_ps(pmat[j4][m-2],
                           pmat[j3][m-1],pmat[j2][m],0);        // mcw -> [ 0, T2,M , T3,M-1,T4,M-2]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
                                                                //makespan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
        _mm_store_ps(res,makespan);
        L[m-2] = res[3];
        // m - 2
        mc = makespan;
        mc = _mm_and_ps(mc,mask0110);                           // mc -> [0 , C2,M, C3,M-1, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , C2,M, C3,M-1]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
        mcw = _mm_set_ps(pmat[j4][m-1],pmat[j3][m],0,0);        // mcw -> [ 0, 0 , T3,M,T4,M-1]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M-1]
        _mm_store_ps(res,makespan);
        L[m-1] = res[3];
        //m - 1
        mc = makespan;
        mc = _mm_and_ps(mc,mask0010);                           // mc -> [0 , 0, C3,M, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , 0, C3,M]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , C3,M , max(C3,M, C4,M-1) ]
        mcw = _mm_set_ps(pmat[j4][m],0,0,0);                    // mcw -> [ 0, 0 , 0,T4,M]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M]
        _mm_store_ps(res,makespan);
        L[m] = res[3];
        previousMachineEndTime[j] = res[0];
        previousMachineEndTime[j+1] = res[1];
        previousMachineEndTime[j+2] = res[2];
        previousMachineEndTime[j+3] = res[3];
        j+=4;

        /* At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }
    if(r4>0)
    {
        for(int i=0;i<(4-r4);i++)
        {
            sol.pop_back();
            previousMachineEndTime.pop_back();
        }
    }

}

