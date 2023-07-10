# Demmin Fieldwork

Welcome to the Demmin project. This fieldwork took place from 29th May to 3rd July in 2023 in two fields near the small 
village of Schmarsow located in north-eastern part of Western-Pomerania, Germany.  In this time, the fields were covered
with Winterwheat, planted in the autumn (usually September / October) of the previous year that will be harvested around
Juna or July. With this, the winter wheat we looked at, was already in a pretty advanced stage of growing. \
In those two fields, in total 78 GPS points were marked for measuring the stage of growing, height, density and seeding
angle of the plants as well as to collect samples and measuring the Above Ground Biomass (AGB).

### Course project definition
This Git is dedicated to build up a linear model to estimate AGB of winter wheat in Western-Pomerania using additional
Sentinel-2 data, as for Landsat, there is no imagery of high enough quality for the period of the fieldwork.

### Study Area
To gain a better understanding of the location and the fields we are in, we present a field map as follows:\
\
![Overview](qgis/overview.png)
\
The map shows the location of the fields (pink) in the bottom right corner in Western-Pomerania (yellow outline) as well
as the fields themselves in the more detailed part of the map. Every field contains three *ESU*s, an agglomeration of 13
individual *SSU* points where measuring will be done. The points of each ESU are arranged in a way that they match the
pixel size of Landsat (30x30m) and Sentinel (10x10m) respectively. However, since as previously mentioned Landsat data
is unavailable, we use all pixels for Sentinel data.

### Measurements
For the ground truth selection of the AGB, at each point, samples were harvested and weighted back in a laboratory using a
scale. Furthermore, the samples were put into a dry oven for 24 hours to acquire their dry weight. Also weights of the 
bowl used for weighing were taken which will be subtracted from the weights afterwards, receiving the pure AGB and Dryed
Biomass. For each point we sampled, the unique identifier consist of the **field name**, the **ESU** and the **SSU**
together.
\
### Sentinel-2 Data
The Sentinel-2 scene was acquired on the 3rd June 2023 being the most cloudless scene available for the period of the
fieldwork. The scene(s) (Landsat excluded in this analysis) were harvested via the Google Earth Engine (GEE) script
[here](https://code.earthengine.google.com/ea8954cb78b05868eca1926ac2fd3bdc?noload=true).

### Scripts


