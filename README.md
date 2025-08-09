# Wind_in_Boulder

This scipt compares surface wind forecasts for a single location (Boulder, CO) from 4 idfferent model sources:

- HRRR (NOAA, 3km resolution, only available in North America, https://rapidrefresh.noaa.gov/hrrr/)
- Meteomatics (commercial downscaled model, US, 1km, https://www.meteomatics.com/en/weather-model-us/)
- tomorrow.io (commercial downscaled model, description of the model is a bit confucing, but think it's 2km, https://www.tomorrow.io/blog/how-does-the-technology-behind-tomorrow-io-work/)
- google models (GraphCast and GenCast, not tested yet, global models, ~25km, https://deepmind.google/science/weathernext/#access-weathernext)

and adds latest observations from:
- Synoptic network (see stations in the Viewer https://viewer.synopticdata.com/map/data/now/air-temperature?layers=#map=11.44/39.9991/-105.2799)
- Metar data from Iowa mesonet (Boulder airport metar not available from Iowa mesonet unfortunately, https://mesonet.agron.iastate.edu/request/asos/1min.phtml)




To do:
------
- get a better idea of how the tomorrow model works
- make GraphCast forecast work
- retreive latest forecast isntead of "old forecast" 6h ago (Goggle model has some latency)
- include the past (historical) forecast hours if available for the individual models
- save forecasts and observations as a database to compare
- compare to more/other ground truth measurements (Get data from the RMHPA station? Metar data from Boudler airport is not in Iowa mesonet)
- extend the comparison (upper level winds? more spatial variability? other ideas?)