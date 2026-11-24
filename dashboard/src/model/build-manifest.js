const EDUCTOR_BUILD={"codeSHA256": "af5707ac4eb62885d9d208cd313818704ed4a81997e3e942e72796c97566a242", "sources": {"exposure/assemblies.js": "74ddd6abfb9dafeb274ec7a229dcb01f127042d7978735d85e8d17270a2a0cb0", "exposure/collisions.js": "480783928f9f4938e64eeab427f232b1734103ea573e396182276c736d7597a6", "exposure/dispersion.js": "401eb2de5eb540d48c2cfb9817801b29f5e6b8c0b0c35997914085768c737bcb", "exposure/geometry.json": "bc79d875e51bdf344b55f230a278079f95da278221422e7d094cc55dbb3cb009", "exposure/model.js": "c480474d489eeaf9756685742b2f8c9868caa209dfa8ad901798d80a1d992abd", "exposure/network.js": "fdd316d6679443c2975eddcdf7886bea05f56c03438c6faee76c99358af9e162", "exposure/numerics.js": "c366f11a66e0ebd24e2b1ed822ca83e769dfcb4959c79d2bccb958014575f8ca", "exposure/population.js": "c29e6414fa916fbaaf1e363f5a8fdcaccf27f5e5e95b38f1db7b5044e7ddce2c", "exposure/runtime.js": "0407d9bcbd438591bfc83dfc980b473a43b9f03bc948c2b54abaa76427bbbb3c", "exposure/worker.js": "646cd58ce19f90a74b79e8e92ead1dc22ee069f3718c325503a16f23895f25e6", "phases/calibration-candidate.json": "2d55517d65f11761b6e685c24170a72a4b661edb0faf2e711889e968251e8252", "phases/walls.json": "87ab6bd67451bedeb826bef9e6f46f7a5ed92f2b7cc90521c18b363f4c8b9a2b", "upgrades/calibration.json": "062a6020a1844932c81a0f115c5898ebfb1e37b8e76e8d36ee1ff5fdc8bb7d79", "upgrades/continuation.js": "08dd4234d74e7adc89595d514d7d2a359dea1b8abd82719168d3c1d96d375d7b", "upgrades/jet.js": "653429c765d91545f140859069bbc8a8dbc0aff6baa48a28311c40f5207e4dd7", "upgrades/physics.js": "4426523bdb8a5b11a502457c84d50d530d9901decb743be64c87dc4ce74ebf05", "upgrades/spatial.js": "e041ed4a4a3c667cb1e472cb37265321420566626da41234392c7dc0bdb54f25", "upgrades/topology.json": "ac17cdfc828c33a27681e63ecb03a0fc8b16bdc101568103896d40baeb6f80b7", "upgrades/walls.js": "16bc58c14373ac9e18e9da658af6c2567d5f5efd3f5664db79f2c7f381d81f11"}, "buildDate": "2026-11-23T06:08:12.124776+00:00"};const EDUCTOR_WALLS={};const EDUCTOR_CANDIDATE={
  "method": "Family-level log-space MAP calibration with additive log discrepancy and Laplace prediction intervals; NOT a sampled Bayesian posterior.",
  "priors": {
    "logFeederSD": 1.5,
    "logMixingSD": 1.5,
    "logDiscrepancySD": 0.1,
    "observationLogSD": 0.05
  },
  "folds": [
    {
      "id": "british35",
      "family": "british",
      "sourceHash": "cfa61df0b09ba35b0b9be170abb80332416b2869def9b1cf6375ce97fea45a83",
      "training": [
        "british25"
      ],
      "parameters": [
        0.8819595252893219,
        0.027996469716078387,
        -0.027734375000000012,
        0.030761718750000007
      ],
      "raw": [
        0.0058742264882276395,
        0.018793895674413825
      ],
      "prediction": [
        0.00571354696047318,
        0.019381012267633273
      ],
      "targets": [
        0.0057894359934826,
        0.02011617896767593
      ],
      "errorsPct": [
        -1.3108191038790462,
        -3.6546040936699553
      ],
      "intervals": [
        [
          0.004936725271004886,
          0.006612605943714465
        ],
        [
          0.016862691474778202,
          0.02227542602436733
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": -0.7317921863771872
    },
    {
      "id": "british25",
      "family": "british",
      "sourceHash": "fc554dcd6db54b484064673d0265a10308083cb349658c945df0f5dd2127c010",
      "training": [
        "british35"
      ],
      "parameters": [
        0.9798110877893225,
        -0.09016759278392161,
        -0.02646484375000001,
        0.036328125
      ],
      "raw": [
        0.005342882969322652,
        0.013607037586503057
      ],
      "prediction": [
        0.005203339054790178,
        0.014110444297480672
      ],
      "targets": [
        0.005086637371566818,
        0.013839213085114885
      ],
      "errorsPct": [
        2.2942795937390237,
        1.9598745296979025
      ],
      "intervals": [
        [
          0.004542703426654457,
          0.0059600495071377236
        ],
        [
          0.012343645440348966,
          0.016130132644726704
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": -0.7621404778898444
    },
    {
      "id": "cad45-sst-coarse-v3",
      "family": "branched",
      "sourceHash": "927cd55e0c520552deacee73707555e03b992b00ebea1f8af09b843eb88e0c35",
      "training": [
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.444140625,
        0.1401376986857903,
        -0.025439453124999952,
        0.004687500000000001
      ],
      "raw": [
        0.01595128757662176,
        0.018801659691344844
      ],
      "prediction": [
        0.01555061361556474,
        0.018889999355481146
      ],
      "targets": [
        0.0169282678493279,
        0.020094052970119774
      ],
      "errorsPct": [
        -8.138187828932875,
        -5.992089383008381
      ],
      "intervals": [
        [
          0.013992567884499684,
          0.01728214476546983
        ],
        [
          0.0169916189639062,
          0.021000475376011254
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": 0.8075397797885335
    },
    {
      "id": "cad45-p2bar-sst-v1",
      "family": "branched",
      "sourceHash": "9c1786665923443c5c31711d70063c7421223403d506481c75daca88ef9a916e",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.36035156249999983,
        -0.03906152006420971,
        -0.05146484374999996,
        -0.015478515625000005
      ],
      "raw": [
        0.012534779965742222,
        0.014435028042224279
      ],
      "prediction": [
        0.011905998325831064,
        0.014213315551803355
      ],
      "targets": [
        0.011817516366821166,
        0.013923453920772641
      ],
      "errorsPct": [
        0.7487356586898475,
        2.0818227480055462
      ],
      "intervals": [
        [
          0.010714918327135234,
          0.013229479852936168
        ],
        [
          0.012782367967039299,
          0.01580445340769887
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": 0.8060136028137729
    },
    {
      "id": "cad45-p5bar-sst-v1",
      "family": "branched",
      "sourceHash": "7cd0fa1025054b864cbe705720e92083fc102de4d778642bdeb2adeb7cc38c1a",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.5121093750000005,
        0.2739267611857903,
        -0.003027343750000109,
        0.022460937500000052
      ],
      "raw": [
        0.02041822991462932,
        0.025168850105834462
      ],
      "prediction": [
        0.020356510384165912,
        0.025740562639505693
      ],
      "targets": [
        0.023706507914533916,
        0.02900162824184617
      ],
      "errorsPct": [
        -14.131130331153486,
        -11.244422468787862
      ],
      "intervals": [
        [
          0.01830385471738099,
          0.022639357742888996
        ],
        [
          0.023159054282688652,
          0.028609828230058212
        ]
      ],
      "covered": [
        false,
        false
      ],
      "parameterCorrelation": 0.8083168219043618
    },
    {
      "id": "cad45-ds10mm-p3bar-sst-inlet-v2",
      "family": "branched",
      "sourceHash": "0aceedafba7ee3c23b2a8fbe0e97f2623875f68c83d7af08af1b9bf2e72fee8a",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.4107421874999999,
        -0.09814355131420971,
        -0.06162109374999997,
        -0.027099609375000028
      ],
      "raw": [
        0.0053537080980751205,
        0.008204080212798207
      ],
      "prediction": [
        0.005033765586194226,
        0.007984738307088488
      ],
      "targets": [
        0.005253446712289514,
        0.008397425346610617
      ],
      "errorsPct": [
        -4.181657074419021,
        -4.914447256011611
      ],
      "intervals": [
        [
          0.004449665404008581,
          0.005694539628513703
        ],
        [
          0.007124301295611592,
          0.008949094540956406
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": 0.7161137077628039
    },
    {
      "id": "cad45-ds40mm-p3bar-sst-v1",
      "family": "branched",
      "sourceHash": "c9ea0563f02814b0c96aec7fbadbb0a4f18098eb5b5b7b5676ad12234a8f530d",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.2150390625,
        -0.01171777006420971,
        -0.073486328125,
        -0.03364257812499999
      ],
      "raw": [
        0.016812196031755756,
        0.01966256814647884
      ],
      "prediction": [
        0.015621032590582031,
        0.019012072194615787
      ],
      "targets": [
        0.01784467558608473,
        0.020988287002385555
      ],
      "errorsPct": [
        -12.461100706345674,
        -9.41579847628894
      ],
      "intervals": [
        [
          0.014020500403531568,
          0.017404276036721313
        ],
        [
          0.017090845213557516,
          0.02114926936712008
        ]
      ],
      "covered": [
        false,
        true
      ],
      "parameterCorrelation": 0.78930227449065
    },
    {
      "id": "cad45-dm10mm-p3bar-sst-outlet-v2",
      "family": "branched",
      "sourceHash": "7c9e64af7b20af363863347d040fe38ab16b25d8cbc1281962c70411fe324abf",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.43945312499999983,
        -0.4316396450642097,
        -0.06455078124999986,
        -0.03813476562500002
      ],
      "raw": [
        0.006139966080001176,
        0.008990338194724262
      ],
      "prediction": [
        0.005756147628551414,
        0.008653948589670009
      ],
      "targets": [
        0.004243058571968738,
        0.007407365214469467
      ],
      "errorsPct": [
        35.66033866651568,
        16.828971423813144
      ],
      "intervals": [
        [
          0.005152789437627529,
          0.0064301551465944275
        ],
        [
          0.007777265272075454,
          0.009629455029848305
        ]
      ],
      "covered": [
        false,
        false
      ],
      "parameterCorrelation": 0.7711785176714785
    },
    {
      "id": "cad45-dm40mm-p3bar-sst-long-v2",
      "family": "branched",
      "sourceHash": "2979b57a7c29098918dc99c2ea4b7f133e6756c5167b6a6e1e627723e0c5427d",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.20039062500000004,
        -0.06249902006420971,
        -0.06235351562499999,
        -0.02426757812500003
      ],
      "raw": [
        0.01673093941800278,
        0.019581311532725864
      ],
      "prediction": [
        0.015719565544502354,
        0.019111840033150338
      ],
      "targets": [
        0.014912497208186094,
        0.01808412526563588
      ],
      "errorsPct": [
        5.412026738708975,
        5.68296642728614
      ],
      "intervals": [
        [
          0.013885679008252735,
          0.01779565412401099
        ],
        [
          0.016939536043502793,
          0.02156271745074301
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": 0.7564467590788049
    },
    {
      "id": "cad90-angle-only-sst-v1",
      "family": "branched",
      "sourceHash": "8a11b5e4385073447ad1e0453f7a42754943668338d73421f51f99421c13d8a0",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "parameters": [
        0.33867187500000007,
        -0.05566308256420971,
        -0.06064453125000002,
        -0.023730468750000042
      ],
      "raw": [
        0.013374668199193427,
        0.016225040313916514
      ],
      "prediction": [
        0.012587672395083486,
        0.0158445450232701
      ],
      "targets": [
        0.013087960487193136,
        0.016229029400963223
      ],
      "errorsPct": [
        -3.822506131487735,
        -2.3691150480650736
      ],
      "intervals": [
        [
          0.011331262761635442,
          0.0139833926420287
        ],
        [
          0.014273293247391305,
          0.01758876544068181
        ]
      ],
      "covered": [
        true,
        true
      ],
      "parameterCorrelation": 0.8130229020323216
    },
    {
      "id": "cad90-d38p1-dj2-p3g-open-sst-v2",
      "family": "branched",
      "sourceHash": "7efd7c8f9eaa33efecbe55cb85aa2d613df04163843414c5850ad3d1389c75df",
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1"
      ],
      "parameters": [
        0.5384765625000003,
        0.2163095736857903,
        0.018847656250000018,
        0.04262695312500001
      ],
      "raw": [
        0.02226467013013308,
        0.024696987668030114
      ],
      "prediction": [
        0.022688286527205902,
        0.02577250516927548
      ],
      "targets": [
        0.019971019018703552,
        0.022623577870895088
      ],
      "errorsPct": [
        13.606053381440052,
        13.918785597707961
      ],
      "intervals": [
        [
          0.020394876249960357,
          0.025239591514638064
        ],
        [
          0.023225754450788408,
          0.028598512229504092
        ]
      ],
      "covered": [
        false,
        false
      ],
      "parameterCorrelation": 0.8180946066007959
    }
  ],
  "fits": {
    "british": {
      "parameters": [
        1.0419204627893222,
        -0.09407384278392161,
        -0.028124999999999983,
        0.02753906249999997
      ],
      "feederResistance": 2.834655641019057,
      "mixingLoss": 0.9102155469091614,
      "logDiscrepancy": [
        -0.028124999999999983,
        0.02753906249999997
      ],
      "training": [
        "british35",
        "british25"
      ],
      "covariance": [
        [
          0.763910677593919,
          -0.3247924572811161,
          0.04124491360514729,
          -0.04752564184773692
        ],
        [
          -0.32479245728111644,
          0.23910022216532806,
          0.005904300189615462,
          0.03686627308173652
        ],
        [
          0.041244913605147177,
          0.005904300189615533,
          0.008777704892428134,
          0.0013001547948194997
        ],
        [
          -0.04752564184773701,
          0.03686627308173653,
          0.0013001547948194867,
          0.006815632429983213
        ]
      ],
      "parameterCorrelation": -0.7599668610370407
    },
    "branched": {
      "parameters": [
        0.369140625,
        -0.01562402006420971,
        -0.047998046875000006,
        -0.014208984375000002
      ],
      "feederResistance": 1.446491002170128,
      "mixingLoss": 0.9844974017491751,
      "logDiscrepancy": [
        -0.047998046875000006,
        -0.014208984375000002
      ],
      "training": [
        "cad45-sst-coarse-v3",
        "cad45-p2bar-sst-v1",
        "cad45-p5bar-sst-v1",
        "cad45-ds10mm-p3bar-sst-inlet-v2",
        "cad45-ds40mm-p3bar-sst-v1",
        "cad45-dm10mm-p3bar-sst-outlet-v2",
        "cad45-dm40mm-p3bar-sst-long-v2",
        "cad90-angle-only-sst-v1",
        "cad90-d38p1-dj2-p3g-open-sst-v2"
      ],
      "covariance": [
        [
          0.06121176769382574,
          0.0563476734679427,
          0.013709016643888962,
          0.011035228241072802
        ],
        [
          0.05634767346794267,
          0.07964592636821281,
          0.016946255330678377,
          0.013679469662845858
        ],
        [
          0.013709016643888957,
          0.01694625533067838,
          0.004014497367063168,
          0.003019937891541659
        ],
        [
          0.011035228241072798,
          0.013679469662845858,
          0.003019937891541659,
          0.0027060701131147277
        ]
      ],
      "parameterCorrelation": 0.8070059644952623
    }
  },
  "summary": {
    "feederRMS": 13.303836876205887,
    "outletRMS": 8.58930696241025,
    "coverage": 0.6818181818181818,
    "passed": false
  },
  "excluded": [
    {
      "id": "cad45-ds10mm-p3bar-sst-v1",
      "reason": "CFD acceptance failed"
    },
    {
      "id": "cad45-dm10mm-p3bar-sst-v1",
      "reason": "CFD acceptance failed"
    },
    {
      "id": "cad45-dm40mm-p3bar-sst-v1",
      "reason": "CFD acceptance failed"
    },
    {
      "id": "annular-aequal-p3bar-sst-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad90-ball-area00-sst-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad90-ball-area25-sst-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad90-ball-area50-sst-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad90-ball-area75-sst-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-ball00deg-p3bar-sst-mesh2-cpu64-v3-postproc-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-ball30deg-p3bar-sst-mesh2-cpu64-v3-postproc-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-ball60deg-p3bar-sst-mesh2-cpu64-v3-postproc-v1",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-n01-axial-fixedmj-sst-window-v2",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-n04-axial-fixedmj-sst-window-v2",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-n04-inward15-fixedmj-sst-budget-v3",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "cad45-n04-inward30-fixedmj-sst-budget-v4",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "coflow-35-sst-v2",
      "reason": "Unmatched valve, annular, multijet or clean-air topology"
    },
    {
      "id": "air-film-35-sst-v3",
      "reason": "CFD acceptance failed"
    },
    {
      "id": "early-dilution-35-sst-v1",
      "reason": "CFD acceptance failed"
    }
  ],
  "limitations": [
    "Only two accepted British nozzle cases; family uncertainty and extrapolation are weakly identified.",
    "Scalar and particle validation are separate from gas-flow cross-validation.",
    "Model discrepancy is for prediction intervals/comparisons, not a correction to conservation-balanced live flow.",
    "Single-mesh gas CFD cannot establish vacuum, powder-loaded or particle accuracy."
  ]
}
;
'use strict';
