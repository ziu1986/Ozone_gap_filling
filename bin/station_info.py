'''
Container module for ionformation about ozone stations.
'''

class station_coord:
    def __init__(self, lat, lon, alt):
        self.alt = alt
        self.lat = lat
        self.lon = lon

    def __repr__(self):
        return("(lat:%3.2f;lon:%3.2f;alt:%3.2f)" % (self.lat, self.lon, self.alt))

    def __str__(self):
        return("Station location (lat:%3.2f;lon:%3.2f;alt:%3.2f)" % (self.lat, self.lon, self.alt))

station_location = {"Barrow":station_coord(71.32,-156.61,11),"Jergul":station_coord(69.45,24.6,255),"Karasjok":station_coord(69.47,25.22,333),"Svanvik":station_coord(69.45,30.03,30), "Esrange":station_coord(67.88,21.07,475),"Janiskoski":station_coord(68.93,28.85,118), "Pallas":station_coord(67.97,24.12,565), "Prestebakke":station_coord(59,11.53,160)}
