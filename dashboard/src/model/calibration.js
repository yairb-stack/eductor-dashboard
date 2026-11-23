const EDUCTOR_CALIBRATION={
  "version": 2,
  "scope": "Historical British body only (revC); 3.5 mm calibration, 2.5 mm holdout. Not transferred to 20Nov CAD.",
  "trainingCase": "british-p3bar-sst-v7",
  "trainingHash": "cfa61df0b09ba35b0b9be170abb80332416b2869def9b1cf6375ce97fea45a83",
  "dataHash": "fa8b102b76bcc0d85cecaaae429f33e4e26e3fee00f33a91bfbaf5bf48abb012",
  "coefficients": {
    "feederResistance": 4.197021050729741,
    "mixingLoss": 0.7413255389701133,
    "turbulenceScale": 0.72334973667736,
    "spreadScale": 0.9497442890891235,
    "developmentD": 1.3397478822494273
  },
  "targets": {
    "flows": {
      "motive_g_s": 6.643176784387476,
      "feeder_g_s": 5.7894359934826,
      "bleed_g_s": 7.680958758302698,
      "outlet_g_s": 20.116178967675932,
      "inlet_sum_g_s": 20.113571536172774,
      "imbalance_pct": 0.012963543040919533
    },
    "peak": {
      "G_s_inv": 3829961.607121514,
      "position_m": [
        0.0071295211617202125,
        0.0018790282577896128,
        -0.0006263191045274916
      ],
      "x_D": 2.037006046205775
    },
    "shearRmsRadius": 0.00265459375756306
  },
  "flowFit": {
    "best": [
      1.4343750000000002,
      -0.2993154263142097
    ],
    "score": 1.1172669398056097e-7,
    "trace": [
      {
        "round": 0,
        "score": 0.0034573262527046586,
        "parameters": [
          1.2,
          -0.2231435513142097
        ]
      },
      {
        "round": 1,
        "score": 0.001174923232386583,
        "parameters": [
          1.5,
          -0.2231435513142097
        ]
      },
      {
        "round": 2,
        "score": 0.0004519762484833406,
        "parameters": [
          1.5,
          -0.3481435513142097
        ]
      },
      {
        "round": 3,
        "score": 0.000020064609015007302,
        "parameters": [
          1.425,
          -0.2856435513142097
        ]
      },
      {
        "round": 4,
        "score": 0.000020064609015007302,
        "parameters": [
          1.425,
          -0.2856435513142097
        ]
      },
      {
        "round": 5,
        "score": 0.000006412256934044864,
        "parameters": [
          1.44375,
          -0.3012685513142097
        ]
      },
      {
        "round": 6,
        "score": 9.639152391928985e-7,
        "parameters": [
          1.4343750000000002,
          -0.3012685513142097
        ]
      },
      {
        "round": 7,
        "score": 1.4488517963372385e-7,
        "parameters": [
          1.4343750000000002,
          -0.2973623013142097
        ]
      },
      {
        "round": 8,
        "score": 1.1172669398056097e-7,
        "parameters": [
          1.4343750000000002,
          -0.2993154263142097
        ]
      }
    ]
  },
  "fieldFit": {
    "best": [
      -0.32386244393873237,
      -0.051562500000000046,
      0.29248144868579035
    ],
    "score": 0.0008133029425385216,
    "trace": [
      {
        "round": 0,
        "score": 0.010573628078886423,
        "parameters": [
          -0.35667494393873245,
          -0.3,
          0.37685644868579027
        ]
      },
      {
        "round": 1,
        "score": 0.003755152818286751,
        "parameters": [
          -0.35667494393873245,
          -0.15,
          0.37685644868579027
        ]
      },
      {
        "round": 2,
        "score": 0.0019397353379613362,
        "parameters": [
          -0.35667494393873245,
          -0.15,
          0.30185644868579026
        ]
      },
      {
        "round": 3,
        "score": 0.001808047517079651,
        "parameters": [
          -0.35667494393873245,
          -0.15,
          0.2643564486857903
        ]
      },
      {
        "round": 4,
        "score": 0.0015586466371131848,
        "parameters": [
          -0.35667494393873245,
          -0.13125,
          0.2643564486857903
        ]
      },
      {
        "round": 5,
        "score": 0.0008775483432170918,
        "parameters": [
          -0.3285499439387324,
          -0.05625000000000004,
          0.29248144868579035
        ]
      },
      {
        "round": 6,
        "score": 0.0008133029425385216,
        "parameters": [
          -0.32386244393873237,
          -0.051562500000000046,
          0.29248144868579035
        ]
      }
    ]
  },
  "limitations": "A fit is not independent validation. Losses are effective coefficients, not measured component K values. The finite mixing-length closure and local scalar distribution are not CFD-resolved."
}
;
'use strict';
