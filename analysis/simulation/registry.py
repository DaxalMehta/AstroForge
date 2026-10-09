SIMULATIONS = {
    ### LowRes simulations
    "N1": {
        "base": "/home/lewis/data/SEEDZ/production/N1/",
        "base2": "/home/jregan/data/Laureate/ProductionsRuns/Catalogues/CalibrationRuns/Normal1/",
        "feedback": True,
        "array_length": 1
    },

    "N2": {
        "base": "/home/lewis/data/SEEDZ/production/N2/",
        "base2": "/home/jregan/data/Laureate/ProductionsRuns/Catalogues/CalibrationRuns/Normal2/",
        "feedback": True,
        "array_length": 1
    },

    "Rarepeak": {
        "base": "/home/lewis/data/SEEDZ/production/Rarepeak_thermal/",
        "base2": "/home/jregan/data/Laureate/ProductionsRuns/Catalogues/CalibrationRuns/Rarepeak_thermal/",
        "feedback": True,
        "array_length": 350
    },

    ### Renaissance Simulations
    ##  "#d62728"-Red "#1f77b4"-Blue "#2ca02c"-Green

    "Rarepeak_NoFeedback": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Rarepeak/NoFeedback/",
        "base2": "/home/jregan/data/Laureate/ProductionsRuns/Catalogues/HighRes/SmallBoxes/Rarepeak/NoFeedback/",
        "feedback": False,
        "array_length": 1 ,
        "import_dir": "/project/home/p200997/daxal/NoFeedback_Elliptical/output/",
        "label": "RP-NoFB",
        "color": "#1f77b4",
        "linestyle": "-",
        "marker": "s"
    },

    "Rarepeak_WeakFeedback": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Rarepeak/WeakFeedback/",
        "base2": "/home/jregan/data/Laureate/ProductionsRuns/Catalogues/HighRes/SmallBoxes/Rarepeak/WeakFeedback/",
        "feedback": True,
        "array_length": 1,
        "import_dir": "/project/home/p200997/daxal/WeakFB_Elliptical/output/",
        "label": "RP-WeakFB",
        "color": "#2ca02c",
        "linestyle": "--",
        "marker": "s"
    },
    "Rarepeak_Feedback": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Rarepeak/Feedback/",
        "base2": None,
        "feedback": True,
        "array_length": 1,
        "import_dir": "/project/home/p200889/daxal/Rarepeak/Feedback/output/",
        "label": "RP-FidFB",
        "color": "#d62728",
        "linestyle": ":",
        "marker": "s"
    },

    "Normal1_NoFeedback": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Normal1/0.5Mpc/NoFeedback/",
        "base2": "/home/jregan/data/Laureate/ProductionsRuns/Catalogues/HighRes/SmallBoxes/Normal1/NoFeedback/",
        "feedback": False,
        "array_length": 1,
        "import_dir": "/project/home/p200997/daxal/Normal1/NoFeedback/output/",
        "label": "N1-NoFB",
        "color": "#1f77b4",
        "linestyle": "-",
        "marker": "o"
    },

    "Normal1_WeakFeedback": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Normal1/0.5Mpc/WeakFeedback/",
        "base2": None,
        "feedback": True,
        "array_length": 1,
        "import_dir": "/project/home/p200889/daxal/Normal1/WeakFeedback/output/",
        "label": "N1-WeakFB",
        "color": "#2ca02c",
        "linestyle": "--",
        "marker": "o"
    },

    "Normal1_Feedback": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Normal1/0.5Mpc/Feedback/",
        "base2": "/home/jregan/data/Laureate/ProductionsRuns/Catalogues/HighRes/SmallBoxes/Normal1/FullFeedback/",
        "feedback": True,
        "array_length": 1,
        "import_dir": "/project/home/p200997/daxal/Normal1/Feedback/output/",
        "label": "N1-FidFB",
        "color": "#d62728",
        "linestyle": ":",
        "marker": "o"
    },

    "Normal2_NoFeedback": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Normal2/0.5Mpc/NoFeedback/",
        "base2": "/home/jregan/data/Laureate/ProductionsRuns/Catalogues/CalibrationRuns/Normal1/NoFeedback/",
        "feedback": False,
        "array_length": 1,
        "import_dir": "/project/home/p200997/daxal/Normal2_thermal/NoFeedback/output/",
        "label": "N2-NoFB",
        "color": "#1f77b4",
        "linestyle": "-",
        "marker": "^"
    },

    "Normal2_WeakFeedback": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Normal2/0.5Mpc/WeakFeedback/",
        "base2": None,
        "feedback": True,
        "array_length": 1,
        "import_dir": "/project/home/p200997/daxal/Normal2_thermal/WeakFeedback/output/",
        "label": "N2-WeakFB",
        "color": "#2ca02c",
        "linestyle": "--",
        "marker": "^"
    },

    "Normal2_Feedback": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Normal2/0.5Mpc/Feedback/",
        "base2": None,
        "feedback": True,
        "array_length": 1,
        "import_dir": "/project/home/p200997/daxal/Normal2_thermal/Feedback/output/",
        "label": "N2-FidFB",
        "color": "#d62728",
        "linestyle": ":",
        "marker": "^"
    },

    "Normal2_NoFeedback_SNe": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Normal2_SNe/0.5Mpc/NoFeedback/",
        "feedback": False,
        "array_length": 1,
        "import_dir": "/project/home/p200997/daxal/Normal2/NoFeedback/output/",
        "label": "N2-NoFB",
        "color": "#2ca02c",
        "linestyle": "-",
        "marker": "o"
    },

    "Normal2_WeakFeedback_SNe": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Normal2_SNe/0.5Mpc/WeakFeedback/",
        "feedback": True,
        "array_length": 1,
        "import_dir": "/project/home/p200997/daxal/Normal2/WeakFeedback/output/",
        "label": "N2-WeakFB",
        "color": "#2ca02c",
        "linestyle": "--",
        "marker": "s"
    },

    "Normal2_Feedback_SNe": {
        "base": "/home/daxal/data/ProductionRuns/Renaissance/Normal2_SNe/0.5Mpc/Feedback/",
        "feedback": True,
        "array_length": 1,
        "import_dir": "/project/home/p200997/daxal/Normal2/Feedback/output/",
        "label": "N2-FidFB",
        "color": "#2ca02c",
        "linestyle": ":",
        "marker": "^"
    },

    ### Nature Paper Simulations Light Seeds
    "L13": {
        "base": "/home/daxal/data/ProductionRuns/NaturePaperData/Level13/",
        "feedback": False,
        "array_length": 350,
        "label": "L13",
        "color": "#A00000"
    },

    "L14": {
        "base": "/home/daxal/data/ProductionRuns/NaturePaperData/Level14/",
        "feedback": False,
        "array_length": 350,
        "label": "L14",
        "color": "#0072B2"
    },

    "L15": {
        "base": "/home/daxal/data/ProductionRuns/NaturePaperData/Level15/",
        "feedback": False,
        "array_length": 350,
        "label": "L15",
        "color": "#3C0550"
    },

    "L15_Feedback": {
        "base": "/home/daxal/data/ProductionRuns/NaturePaperData/Level15_feedback/",
        "feedback": True,
        "array_length": 350,
        "label": "L15_BHFB",
        "color": "#E69F00"
    },

    ### Isolated Galaxy - Light Seeds

    "Feedback": {
        "base": "/home/daxal/data/ProductionRuns/IsolatedGalaxy/Feedback/",
        "feedback": True,
        "array_length": 1,
        "label": "SN_FB"
    },
    
    "NoFeedback": {
        "base": "/home/daxal/data/ProductionRuns/IsolatedGalaxy/NoFeedback/",
        "feedback": False,
        "array_length": 1,
        "label": "NoSN_FB"
    }

}