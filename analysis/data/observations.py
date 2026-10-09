import numpy as np

class ObservationData:
    def __init__(self):
        self.Mailino()
        self.Juodzbalis()
        self.SpecialLRDs()
        self.Taylor()

    def Mailino(self):
        self.maiolino = {}
        self.maiolino[10013704] = {
            "MBH": [7.5, 0.3, 0.3],
            "MStar": [8.88, 0.66, 0.66],
            "f_edd": 0.06,
            "Redshift": 5.9193,
            "Lbol" : 44.29           
        }
        self.maiolino[8084] = {
            "MBH": [7.25, 0.3, 0.3],
            "MStar": [8.45, 0.03, 0.03],
            "f_edd": 0.16,
            "Redshift": 4.6482,
            "Lbol" : 44.25           
        }
        self.maiolino[1093] = {
            "MBH": [7.36, 0.3, 0.3],
            "MStar": [8.34, 0.2, 0.2],
            "f_edd": 0.2,
            "Redshift": 5.5951,
            "Lbol" : 44.32           
        }
        self.maiolino[3608] = {
            "MBH": [6.82, 0.3, 0.3],
            "MStar": [8.38, 0.11, 0.15],
            "f_edd": 0.11,
            "Redshift": 5.26894,
            "Lbol" : None           
        }
        self.maiolino[11836] = {
            "MBH": [7.12, 0.3, 0.3],
            "MStar": [7.79, 0.3, 0.3],
            "f_edd": 0.2,
            "Redshift": 4.40935,
            "Lbol" : 44.11           
        }
        self.maiolino[20621] = {
            "MBH": [7.3, 0.3, 0.3],
            "MStar": [8.06, 0.7, 0.7],
            "f_edd": 0.18,
            "Redshift": 4.68123,
            "Lbol" : 44.17           
        }
        self.maiolino[73488] = {
            "MBH": [7.71, 0.3, 0.3],
            "MStar": [9.78, 0.2, 0.2],
            "f_edd": 0.16,
            "Redshift": 4.1332,
            "Lbol" : 45.22          
        }
        self.maiolino[77652] = {
            "MBH": [6.86, 0.3, 0.3],
            "MStar": [7.87, 0.16, 0.28],
            "f_edd": 0.38,
            "Redshift": 5.22943,
            "Lbol" : 44.11           
        }
        self.maiolino[61888] = {
            "MBH": [7.22, 0.3, 0.3],
            "MStar": [8.11, 0.92, 0.92],
            "f_edd": 0.32,
            "Redshift": 5.87461,
            "Lbol" : 44.38           
        }
        self.maiolino[62309] = {
            "MBH": [6.56, 0.3, 0.3],
            "MStar": [8.12, 0.12, 0.13],
            "f_edd": 0.39,
            "Redshift": 5.17241,
            "Lbol" : 43.56           
        }
        self.maiolino[53757] = {
            "MBH": [7.69, 0.3, 0.3],
            "MStar": [10.18, 0.13, 0.12],
            "f_edd": 0.05,
            "Redshift": 4.4480,
            "Lbol" : 44.29           
        }
        self.maiolino[954] = {
            "MBH": [7.9, 0.3, 0.3],
            "MStar": [10.66, 0.09, 0.1],
            "f_edd": 0.42,
            "Redshift": 6.76026,
            "Lbol" : 45.17           
        }
    
    def Juodzbalis(self):
        self.juodzbalis = {}
        self.juodzbalis[30148179] = {
            "MBH": [7.12, 0.3, 0.3],
            "MStar": [8.95, 0.58, 0.58],
            "f_edd": 0.11,
            "Redshift": 5.922,
            "Lbol" : 44.25           
        }
        self.juodzbalis[210600] = {
            "MBH": [7.41, 0.3, 0.3],
            "MStar": [8.40, 0.61, 0.61],
            "f_edd": 0.57,
            "Redshift": 6.306,
            "Lbol" : 44.29           
        }
        self.juodzbalis[209777] = {
            "MBH": [8.90, 0.3, 0.3],
            "MStar": [None, None, None],
            "f_edd": 0.27,
            "Redshift": 3.709,
            "Lbol" : 45.42           
        }
        self.juodzbalis[179198] = {
            "MBH": [7.23, 0.3, 0.3],
            "MStar": [8.44, 0.05, 0.05],
            "f_edd": 0.04,
            "Redshift": 3.830,
            "Lbol" : 43.92           
        }
        self.juodzbalis[172975] = {
            "MBH": [7.25, 0.3, 0.3],
            "MStar": [8.98, 0.14, 0.14],
            "f_edd": 0.15,
            "Redshift": 4.741,
            "Lbol" : 44.07           
        }
        self.juodzbalis[159717] = {
            "MBH": [7.44, 0.3, 0.3],
            "MStar": [None, None, None],
            "f_edd": 0.38,
            "Redshift": 5.077,
            "Lbol" : 45.13           
        }
        self.juodzbalis[159438] = {
            "MBH": [6.47, 0.3, 0.3],
            "MStar": [8.35, 0.13, 0.13],
            "f_edd": 0.35,
            "Redshift": 2.239,
            "Lbol" : 44.11           
        }
        self.juodzbalis[49729] = {
            "MBH": [7.67, 0.3, 0.3],
            "MStar": [None, None, None],
            "f_edd": 0.115,
            "Redshift": 3.189,
            "Lbol" : 44.83           
        }
        self.juodzbalis[38562] = {
            "MBH": [7.51, 0.3, 0.3],
            "MStar": [9.76, 0.09, 0.09],
            "f_edd": 0.12,
            "Redshift": 4.822,
            "Lbol" : 44.70           
        }
        self.juodzbalis[29648] = {
            "MBH": [6.81, 0.3, 0.3],
            "MStar": [9.71, 0.01, 0.01],
            "f_edd": 0.11,
            "Redshift": 2.960,
            "Lbol" : 43.90           
        }
        self.juodzbalis[17341] = {
            "MBH": [6.76, 0.3, 0.3],
            "MStar": [8.54, 0.03, 0.03],
            "f_edd": 0.15,
            "Redshift": 3.598,
            "Lbol" : 44.01           
        }
        self.juodzbalis[13329] = {
            "MBH": [6.86, 0.3, 0.3],
            "MStar": [9.52, 0.14, 0.14],
            "f_edd": 0.14,
            "Redshift": 3.936,
            "Lbol" : 44.11           
        }
        self.juodzbalis[9598] = {
            "MBH": [6.48, 0.3, 0.3],
            "MStar": [9.08, 0.38, 0.38],
            "f_edd": 0.18,
            "Redshift": 3.324,
            "Lbol" : 43.85           
        }
        self.juodzbalis[2916] = {
            "MBH": [7.05, 0.3, 0.3],
            "MStar": [8.98, 0.97, 0.97],
            "f_edd": 0.06,
            "Redshift": 3.664,
            "Lbol" : 43.91           
        }
        self.juodzbalis[1001830] = {
            "MBH": [8.61, 0.37, 0.38],
            "MStar": [8.92, 0.31, 0.3],
            "Redshift": 6.68,
            "f_edd": 0.024,
            "Lbol": None

        }

    def SpecialLRDs(self):
        self.speciallrds = {}
        self.speciallrds["Maiolino 2024"] = {
            "MBH": [6.2, 0.3, 0.3],
            "MStar": [8.90, 0.6, 0.3],
            "f_edd": 5.5,
            "Redshift": 10.6034,
            "Lbol" : 45,
            "color": "magenta",
        }
        self.speciallrds["Kovacs 2024"] = {
            "MBH": [7.90, 0.40, 0.34],
            "MStar": [7.69, 0.28, 0.39],
            "color": "red",
            "f_edd": None,
            "Redshift": 10,
            "Lbol": 46,
        }
        self.speciallrds["Bogdan 2024"] = {
            "MBH": [7.60, 1.00, 0.18],
            "MStar": [7.60, 0.30, 0.30],
            "Redshift": 10.3,
            "Lbol": 45.70,
            "f_edd": None,
            "color": "bisque",
        }
        self.speciallrds["Taylor 2025"] = {
            "MBH": [7.58, 0.15, 0.15],
            "MStar": [8.9, 0.1, 0.1],
            "Redshift": 9.288,
            "Lbol": None,
            "f_edd": None,
            "color": "darkred",
        }
        self.speciallrds["Ortiz 2025"] = {
            "MBH": [7.2, 0.04, 0.04],
            "MStar": [8.5, 0.27, 0.14],
            "Redshift": 12.34,
            "color": "cyan",
        }

    def Taylor(self):
        self.taylor = {}
        self.taylor[0] = {
            "MBH" : 6.25,
            "Phi" : [258, 125, 113]
        }
        self.taylor[1] = {
            "MBH" : 6.75,
            "Phi" : [276, 83.7, 59.8]
        }
        self.taylor[2] = {
            "MBH" : 7.25,
            "Phi" : [113, 32.5, 24.3]
        }
        self.taylor[3] = {
            "MBH" : 7.75,
            "Phi" : [36.1, 14.8, 10.5]
        }
        self.taylor[4] = {
            "MBH" : 8.25,
            "Phi" : [7.67, 7.51, 4.17]
        }

    def extract_data(self, dataset, field):
        data = np.array([dataset[k][field] for k in dataset], dtype=float)
        out = {
            "val": data[:, 0],
            "err": np.stack((data[:, 1], data[:, 2]))
        }
        return out

    def Reines(self, mstellar):
        alpha, beta = 7.45, 1.05
        alpha_err, beta_err = 0.08, 0.11

        mbh = alpha + beta * np.log10(mstellar/1e11)
        #mbh_err = mbh * np.sqrt((alpha_err / alpha)**2 + (beta_err / beta)**2 + (0.1 / np.log(10))**2)
        mbh_err = 0.55

        return mbh, mbh_err

    def Pacucci(self, mstellar):
        alpha, beta = -2.43, 1.06
        alpha_err, beta_err = 0.83, 0.09
        mbh = alpha + beta * np.log10(mstellar)
        #mbh_err = mbh * np.sqrt((alpha_err / alpha)**2 + (beta_err / beta)**2 + (0.1 / np.log(10))**2)
        mbh_err = 0.69
        return mbh, mbh_err

    def Li(self, mstellar):
        alpha, beta = -3.75, 1.07
        mbh = alpha + beta * np.log10(mstellar)
        mbh_err = 0.91

        return mbh, mbh_err

    def lisa_events(self):

        z = np.array([10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0])

        curves = {
            "BACH": {
                "data": np.array([132, 170, 201, 94.2, 135, 112, 38.8, 201, 249, 201, np.nan]),
                "color": "#7B2CBF"
            },

            "CAT": {
                "data": np.array([1.27, 1.27, 2.25, 1.71, 0.871, 1.98, 1.44, 2.20, 2.02, 1.27, 5.1]),
                "color": "#E76F51"
            },

            "Renaissance": {
                "data": np.array([np.nan, np.nan, np.nan, 14.4, 19.8, 16.7, 12.4, 8.86, np.nan, np.nan, np.nan]),
                "color": "#E9C46A"
            },

            "Romulus25": {
                "data": np.array([0.303, 0.297, 0.256, 0.128, 0.154, 0.168, 0.0786, 0.0950, np.nan, np.nan, np.nan]),
                "color": "#2A9D8F"
            },
            "L-Galaxies": {
                "data": np.array([0.536, 0.330, 0.199, 0.110, 0.0664, 0.0435, 0.0222, 0.0113, 6.38e-3, 1.69e-3, 1.7e-3]),
                "color": "#A23B72"
            }
        }

        return z, curves