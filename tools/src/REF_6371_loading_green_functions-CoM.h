/*
    Love Numbers and Green Functions by Pascal Gegout
    Geosciences Environnement Toulouse - April 2016 -
 
    Notations:
    - PsiDeg: Angular Distance from load in degrees,
    - PsiRad: Angular Distance from load in radians,
 
    Green's functions (following Farrell's units):
    - Gu: Vertical   Displacement (1e12*a*PsiRad meters),
    - Gv: Horizontal Displacement (1e12*a*PsiRad meters),
    - Gf: Gravitational Potential (1e12*a*PsiRad meters),
    - Gg: Gravity Perturbation (1e18.a.PsiRad m.s**(-2)),
    - Gs: Self-Attraction plus Loading (Gs=Gf-Gu meters),
    - Gt: Gravity Tilt  (1e12*a*a*PsiRad*PsiRad radians),
    - Gh: Horizontal Gradient of Self-Attraction+Loading 
          Gh = -a * Gt
 
    Constants:
    - Earth's Radius  a = 6371000e00 m
    - Earth's Mass    M = 5.97219e24 kg
    - Earth's Gravity g = 9.81e0 m.s**(-2)
 
    Remarks:
    - Reference: Center of Mass (all masses, load included).
    - Compressibility included in all Loading Love Numbers.
 
    References:
    - Loading Love Numbers by P. Gegout using Earth's model
      REF_6371, DATA included in the REAR software package,
      (D. Melini, P. Gegout, M. King, B.Marzeion, G. Spada,
       EOS 2015 and https://github.com/danielemelini/rear).
*/
 
/*
        Angular      Radial       Horizontal   Gravitational     Vertical     Loading+Self       Gravity      Horizontal
       Distance   Displacement   Displacement    Potential        Gravity     -Attraction         Tilt       Gradient of
        PsiDeg         Gu             Gv             Gf             Gg           LSA  Gs           Gt          LSA  Gh  
        degree         m              m              m           m.s**(-2)           m           radians          m     
          Norm:     x 1e12         x 1e12         x 1e12         x 1e18         x 1e12         x 1e12         x 1e12    
          Norm:     x a.PsiRad     x.a.PsiRad     x a.PsiRad     x a.PsiRad     x a.PsiRad     x a.PsiRad     x - PsiRad
          Norm:                                                                                x a.PsiRad     x a.PsiRad
*/
  {   0.0000001,    -42.238780,    -12.858080,      6.794509,   -103.201968,     49.033289,     49.033945,     49.033945},
  {   0.0000005,    -42.238695,    -12.858080,      6.794508,   -103.201781,     49.033204,     49.033944,     49.033944},
  {   0.0000010,    -42.238589,    -12.858080,      6.794507,   -103.201549,     49.033097,     49.033944,     49.033944},
  {   0.0000020,    -42.238377,    -12.858080,      6.794505,   -103.201083,     49.032882,     49.033944,     49.033944},
  {   0.0000040,    -42.237953,    -12.858080,      6.794501,   -103.200152,     49.032454,     49.033943,     49.033943},
  {   0.0000060,    -42.237529,    -12.858080,      6.794497,   -103.199221,     49.032025,     49.033942,     49.033942},
  {   0.0000080,    -42.237104,    -12.858081,      6.794492,   -103.198290,     49.031597,     49.033942,     49.033942},
  {   0.0000100,    -42.236680,    -12.858081,      6.794488,   -103.197359,     49.031168,     49.033941,     49.033941},
  {   0.0000200,    -42.234559,    -12.858082,      6.794467,   -103.192703,     49.029026,     49.033937,     49.033937},
  {   0.0000400,    -42.230316,    -12.858083,      6.794425,   -103.183392,     49.024741,     49.033930,     49.033930},
  {   0.0000600,    -42.226073,    -12.858085,      6.794383,   -103.174081,     49.020456,     49.033923,     49.033923},
  {   0.0000800,    -42.221830,    -12.858087,      6.794341,   -103.164770,     49.016171,     49.033916,     49.033916},
  {   0.0001000,    -42.217587,    -12.858088,      6.794299,   -103.155459,     49.011886,     49.033908,     49.033908},
  {   0.0002000,    -42.196372,    -12.858088,      6.794088,   -103.108904,     48.990460,     49.033872,     49.033872},
  {   0.0004000,    -42.153944,    -12.858067,      6.793667,   -103.015796,     48.947611,     49.033798,     49.033798},
  {   0.0006000,    -42.111516,    -12.858019,      6.793247,   -102.922692,     48.904763,     49.033721,     49.033721},
  {   0.0008000,    -42.069092,    -12.857941,      6.792827,   -102.829594,     48.861919,     49.033639,     49.033639},
  {   0.0010000,    -42.026670,    -12.857835,      6.792409,   -102.736504,     48.819078,     49.033552,     49.033552},
  {   0.0020000,    -41.814628,    -12.856874,      6.790337,   -102.271228,     48.604964,     49.032983,     49.032983},
  {   0.0040000,    -41.391128,    -12.852794,      6.786374,   -101.342134,     48.177502,     49.030790,     49.030790},
  {   0.0060000,    -40.968765,    -12.845811,      6.782739,   -100.415608,     47.751504,     49.026886,     49.026886},
  {   0.0080000,    -40.547709,    -12.835895,      6.779428,    -99.491578,     47.327137,     49.021564,     49.021564},
  {   0.0100000,    -40.127939,    -12.823040,      6.776351,    -98.569628,     46.904290,     49.015202,     49.015202},
  {   0.0200000,    -38.050222,    -12.715559,      6.762567,    -94.000867,     44.812789,     48.950339,     48.950339},
  {   0.0300000,    -36.021810,    -12.539174,      6.750634,    -89.532139,     42.772443,     48.798470,     48.798470},
  {   0.0400000,    -34.061394,    -12.299745,      6.739867,    -85.203291,     40.801261,     48.524662,     48.524662},
  {   0.0500000,    -32.186021,    -12.004947,      6.729932,    -81.050382,     38.915953,     48.103144,     48.103144},
  {   0.0600000,    -30.410049,    -11.663799,      6.720624,    -77.103808,     37.130674,     47.518760,     47.518760},
  {   0.0700000,    -28.744670,    -11.286145,      6.711808,    -73.387392,     35.456478,     46.767328,     46.767328},
  {   0.0800000,    -27.197729,    -10.882135,      6.703382,    -69.918038,     33.901111,     45.854997,     45.854997},
  {   0.0900000,    -25.773768,    -10.461754,      6.695271,    -66.705820,     32.469040,     44.796790,     44.796790},
  {   0.1000000,    -24.474246,    -10.037798,      6.687418,    -63.754439,     31.161663,     43.614606,     43.614606},
  {   0.1500000,    -19.727497,     -8.051574,      6.650686,    -52.693261,     26.378183,     36.813983,     36.813983},
  {   0.2000000,    -17.251659,     -6.665981,      6.616332,    -46.497392,     23.867991,     30.711052,     30.711052},
  {   0.2500000,    -16.088106,     -5.928556,      6.583290,    -43.206738,     22.671396,     26.648410,     26.648410},
  {   0.3000000,    -15.558929,     -5.625956,      6.551268,    -41.404569,     22.110196,     24.426336,     24.426336},
  {   0.3500000,    -15.297313,     -5.574536,      6.520187,    -40.302638,     21.817501,     23.400875,     23.400875},
  {   0.4000000,    -15.131679,     -5.634379,      6.490017,    -39.514592,     21.621696,     23.031665,     23.031665},
  {   0.4500000,    -14.989752,     -5.734099,      6.460722,    -38.865069,     21.450473,     22.978881,     22.978881},
  {   0.5000000,    -14.844658,     -5.835280,      6.432257,    -38.277591,     21.276915,     23.055960,     23.055960},
  {   0.5500000,    -14.688760,     -5.920250,      6.404574,    -37.721087,     21.093335,     23.169470,     23.169470},
  {   0.6000000,    -14.521858,     -5.983864,      6.377622,    -37.183263,     20.899481,     23.276903,     23.276903},
  {   0.6500000,    -14.346164,     -6.025144,      6.351352,    -36.659938,     20.697516,     23.361794,     23.361794},
  {   0.7000000,    -14.164295,     -6.045985,      6.325717,    -36.149827,     20.490012,     23.420008,     23.420008},
  {   0.7500000,    -13.978581,     -6.049368,      6.300673,    -35.652446,     20.279254,     23.452687,     23.452687},
  {   0.8000000,    -13.790890,     -6.038034,      6.276181,    -35.167829,     20.067071,     23.462851,     23.462851},
  {   0.8500000,    -13.602655,     -6.014641,      6.252204,    -34.695528,     19.854859,     23.453839,     23.453839},
  {   0.9000000,    -13.414950,     -5.981789,      6.228708,    -34.235530,     19.643658,     23.428701,     23.428701},
  {   0.9500000,    -13.228572,     -5.941109,      6.205664,    -33.787215,     19.434236,     23.390086,     23.390086},
  {   1.0000000,    -13.044107,     -5.894608,      6.183042,    -33.350423,     19.227149,     23.340173,     23.340173},
  {   1.1000000,    -12.682501,     -5.788686,      6.138969,    -32.509311,     18.821470,     23.213516,     23.213516},
  {   1.2000000,    -12.332274,     -5.672106,      6.096309,    -31.708868,     18.428583,     23.059975,     23.059975},
  {   1.3000000,    -11.994424,     -5.550052,      6.054908,    -30.945950,     18.049332,     22.887160,     22.887160},
  {   1.4000000,    -11.669315,     -5.425826,      6.014638,    -30.217750,     17.683953,     22.700013,     22.700013},
  {   1.5000000,    -11.356973,     -5.301536,      5.975385,    -29.521842,     17.332358,     22.501711,     22.501711},
  {   1.6000000,    -11.057240,     -5.178499,      5.937053,    -28.856136,     16.994293,     22.294702,     22.294702},
  {   1.7000000,    -10.769822,     -5.057520,      5.899556,    -28.218770,     16.669377,     22.081552,     22.081552},
  {   1.8000000,    -10.494307,     -4.939079,      5.862820,    -27.608003,     16.357127,     21.864640,     21.864640},
  {   1.9000000,    -10.230226,     -4.823442,      5.826779,    -27.022221,     16.057004,     21.645215,     21.645215},
  {   2.0000000,     -9.977120,     -4.710720,      5.791375,    -26.460008,     15.768495,     21.423550,     21.423550},
  {   2.2500000,     -9.389658,     -4.441711,      5.705323,    -25.149971,     15.094981,     20.861581,     20.861581},
  {   2.5000000,     -8.862138,     -4.190518,      5.622261,    -23.963947,     14.484399,     20.295347,     20.295347},
  {   2.7500000,     -8.388486,     -3.956239,      5.541664,    -22.887924,     13.930151,     19.737206,     19.737206},
  {   3.0000000,     -7.963450,     -3.737573,      5.463121,    -21.910418,     13.426572,     19.182405,     19.182405},
  {   3.2500000,     -7.582457,     -3.533687,      5.386303,    -21.021723,     12.968760,     18.643440,     18.643440},
  {   3.5000000,     -7.240944,     -3.343554,      5.310945,    -20.212512,     12.551889,     18.115321,     18.115321},
  {   3.7500000,     -6.935522,     -3.166458,      5.236835,    -19.475936,     12.172356,     17.604942,     17.604942},
  {   4.0000000,     -6.662279,     -3.001685,      5.163797,    -18.804264,     11.826076,     17.114525,     17.114525},
  {   4.3333333,     -6.343143,     -2.799886,      5.067841,    -17.999880,     11.410984,     16.487278,     16.487278},
  {   4.6666667,     -6.069280,     -2.617219,      4.973276,    -17.287291,     11.042555,     15.898802,     15.898802},
  {   5.0000000,     -5.834856,     -2.452328,      4.879895,    -16.655010,     10.714751,     15.356553,     15.356553},
  {   5.5000000,     -5.546113,     -2.234974,      4.741712,    -15.835556,     10.287825,     14.601504,     14.601504},
  {   6.0000000,     -5.319239,     -2.049856,      4.605464,    -15.144338,      9.924703,     13.940059,     13.940059},
  {   6.5000000,     -5.142238,     -1.892761,      4.470895,    -14.557855,      9.613133,     13.372612,     13.372612},
  {   7.0000000,     -5.005113,     -1.759592,      4.337838,    -14.056442,      9.342951,     12.879231,     12.879231},
  {   7.5000000,     -4.899392,     -1.646744,      4.206185,    -13.623510,      9.105576,     12.460475,     12.460475},
  {   8.0000000,     -4.818486,     -1.550924,      4.075869,    -13.245964,      8.894355,     12.094423,     12.094423},
  {   8.5000000,     -4.756948,     -1.469251,      3.946851,    -12.913016,      8.703799,     11.788189,     11.788189},
  {   9.0000000,     -4.710675,     -1.399115,      3.819111,    -12.616387,      8.529786,     11.532699,     11.532699},
  {   9.5000000,     -4.676040,     -1.338297,      3.692639,    -12.348848,      8.368679,     11.312892,     11.312892},
  {  10.0000000,     -4.650428,     -1.285029,      3.567433,    -12.105010,      8.217861,     11.148128,     11.148128},
  {  11.0000000,     -4.618511,     -1.194640,      3.320831,    -11.672461,      7.939342,     10.864834,     10.864834},
  {  12.0000000,     -4.602319,     -1.117904,      3.079345,    -11.292380,      7.681663,     10.679341,     10.679341},
  {  13.0000000,     -4.594314,     -1.047968,      2.843023,    -10.948540,      7.437337,     10.544900,     10.544900},
  {  14.0000000,     -4.589989,     -0.980268,      2.611905,    -10.630704,      7.201894,     10.445110,     10.445110},
  {  15.0000000,     -4.586059,     -0.911911,      2.386024,    -10.331304,      6.972082,     10.385548,     10.385548},
  {  17.5000000,     -4.570120,     -0.727223,      1.844366,     -9.641971,      6.414485,     10.255120,     10.255120},
  {  20.0000000,     -4.534070,     -0.515146,      1.335767,     -9.006330,      5.869837,     10.187208,     10.187208},
  {  22.5000000,     -4.472288,     -0.272127,      0.860285,     -8.405936,      5.332573,     10.145856,     10.145856},
  {  25.0000000,     -4.381131,      0.000993,      0.417861,     -7.828976,      4.798992,     10.119206,     10.119206},
  {  27.5000000,     -4.260436,      0.303052,      0.008292,     -7.271335,      4.268728,     10.082562,     10.082562},
  {  30.0000000,     -4.109433,      0.634263,     -0.368774,     -6.728524,      3.740659,     10.037856,     10.037856},
  {  35.0000000,     -3.726832,      1.374906,     -1.027619,     -5.695193,      2.699213,      9.892045,      9.892045},
  {  40.0000000,     -3.254288,      2.206082,     -1.564525,     -4.744784,      1.689763,      9.575008,      9.575008},
  {  50.0000000,     -2.138450,      4.042754,     -2.303228,     -3.181656,     -0.164778,      8.433850,      8.433850},
  {  60.0000000,     -0.916709,      5.924735,     -2.649453,     -2.144516,     -1.732745,      6.752491,      6.752491},
  {  70.0000000,      0.342151,      7.655604,     -2.661103,     -1.614604,     -3.003254,      4.956585,      4.956585},
  {  80.0000000,      1.643453,      9.097257,     -2.383359,     -1.480443,     -4.026813,      3.348583,      3.348583},
  {  90.0000000,      3.020882,     10.168624,     -1.852091,     -1.611109,     -4.872973,      2.141944,      2.141944},
  { 100.0000000,      4.509615,     10.825404,     -1.100520,     -1.885966,     -5.610135,      1.414797,      1.414797},
  { 110.0000000,      6.125994,     11.044325,     -0.165545,     -2.215297,     -6.291539,      1.046673,      1.046673},
  { 120.0000000,      7.859488,     10.812487,      0.908045,     -2.544462,     -6.951443,      0.932662,      0.932662},
  { 130.0000000,      9.677491,     10.122431,      2.065708,     -2.839583,     -7.611783,      0.969716,      0.969716},
  { 140.0000000,     11.517863,      8.972418,      3.243524,     -3.101012,     -8.274339,      1.039884,      1.039884},
  { 150.0000000,     13.309648,      7.366336,      4.370332,     -3.329740,     -8.939316,      0.992483,      0.992483},
  { 160.0000000,     14.969619,      5.317985,      5.371082,     -3.538787,     -9.598536,      0.842983,      0.842983},
  { 170.0000000,     16.412598,      2.850426,      6.171048,     -3.745159,    -10.241550,      0.518020,      0.518020},
  { 172.0000000,     16.668669,      2.309117,      6.300620,     -3.785539,    -10.368048,      0.354720,      0.354720},
  { 174.0000000,     16.911050,      1.753461,      6.418809,     -3.829279,    -10.492241,      0.296019,      0.296019},
  { 176.0000000,     17.135505,      1.182929,      6.525087,     -3.882255,    -10.610418,      0.144209,      0.144209},
  { 177.0000000,     17.251328,      0.892675,      6.573613,     -3.892853,    -10.677715,      0.137998,      0.137998},
  { 178.0000000,     17.355204,      0.598413,      6.618951,     -3.918227,    -10.736253,      0.114823,      0.114823},
  { 179.0000000,     17.454613,      0.300547,      6.661105,     -3.947714,    -10.793508,     -0.010718,     -0.010718},
  { 179.2000000,     17.486697,      0.240680,      6.669165,     -3.930296,    -10.817532,     -0.096814,     -0.096814},
  { 179.4000000,     17.498514,      0.180661,      6.677064,     -3.950179,    -10.821449,     -0.005394,     -0.005394},
  { 179.6000000,     17.509598,      0.120548,      6.684831,     -3.971070,    -10.824768,      0.046302,      0.046302},
  { 179.7000000,     17.527896,      0.090470,      6.688683,     -3.957706,    -10.839214,     -0.001524,     -0.001524},
  { 179.8000000,     17.550026,      0.060363,      6.692508,     -3.937328,    -10.857518,     -0.059338,     -0.059338},
  { 179.9000000,     17.568993,      0.030208,      6.696297,     -3.922866,    -10.872696,     -0.062130,     -0.062130},
  { 180.0000000,     17.581875,     -0.000000,      6.700043,     -3.918126,    -10.881833,     -0.000000,     -0.000000} 