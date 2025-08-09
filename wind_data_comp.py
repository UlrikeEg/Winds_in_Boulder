import requests
import datetime
from datetime import datetime, timedelta, UTC
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import tempfile
import os
import pygrib
import ee   # for Google Earth Engine API
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
from Get_obs_data import *
from access_keys import *       # This includes all access keys for model data and synoptoic observation data, put in your own


#%% Definitions 


# Coordinates for Boulder, CO (paragliding hill)
LAT, LON = 40.056731120064015, -105.30018040686129


#Model initiation time
START = datetime.now(UTC).replace(minute=0, second=0, microsecond=0).replace(tzinfo=None)   
# Find the last synoptic time before or equal to START
synoptic_hours = [0, 6, 12, 18]
def last_synoptic_time(start_time):
    for offset in sorted(synoptic_hours, reverse=True):
        candidate = start_time.replace(hour=offset)
        if candidate <= start_time:
            return candidate
    return (start_time - timedelta(days=1)).replace(hour=18)
START = last_synoptic_time(START) - timedelta(hours = 6)   # !!! to be sure that data is available. Replace by last available data.
HOURS =  45   # maximum forecast time for HRRR and Google



    



#%% Functions

### HRRR

def download_hrrr_grib2(forecast_hour, variable="sfc"):
    """
    Download HRRR GRIB2 file for a given forecast hour.
    Example path: https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod/hrrr.20250623/conus/hrrr.t19z.wrfsfcf02.grib2
                  https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod/hrrr.20250624/conus/hrrr.t00z.wrfsfcf02.grib2
    """
    base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod"
    date = START.strftime("%Y%m%d")
    hour = START.strftime("%H")
    fhr = f"{forecast_hour:02d}"
    url = f"{base_url}/hrrr.{date}/conus/hrrr.t{hour}z.wrfsfcf{fhr}.grib2"
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        fd, tmp_path = tempfile.mkstemp(suffix=".grib2")
        with os.fdopen(fd, "wb") as out:
            out.write(response.content)
        return tmp_path
    else:
        print (response)
        return None

def extract_hrrr_variable(grib_path, variable, level=None):
    """
    Extract a variable from HRRR GRIB2 file at Boulder, CO using xarray.
    variable: 'TMP' for temperature, 'UGRD'/'VGRD' for winds.
    level: 'surface' or height in meters AGL (e.g., 1000)
    """
    
    # import cfgrib
    # grbs = pygrib.open(tmp_path)
    # for grb in grbs:
    #     print(grb)
    # for grb in grbs:
    #     print(grb)  # Print variable info (optional)
    #     lats, lons = grb.latlons()
    #     print("Latitude shape:", lats.shape)
    #     print("Longitude shape:", lons.shape)
    #     print("Latitude range: {:.2f} to {:.2f}".format(lats.min(), lats.max()))
    #     print("Longitude range: {:.2f} to {:.2f}".format(lons.min(), lons.max()))
    #     break  # Remove this if you want to check all variables
    # grbs.close()

    
    
    grbs = pygrib.open(grib_path)
    if level == 'surface':
        msgs = grbs.select(name=variable, typeOfLevel='surface')
    elif isinstance(level, int):
        if level <100:  # use meters
            # HRRR gives wind at height above ground (hagl)
            msgs = grbs.select(name=variable, typeOfLevel='heightAboveGround', level=level)
        else:   # use hPa
            msgs = grbs.select(name=variable, typeOfLevel='isobaricInhPa', level=level)   
    else:
        msgs = grbs.select(name=variable)
    # Find nearest point
    lats, lons = msgs[0].latlons()
    idx = np.argmin((lats-LAT)**2 + (lons-LON)**2)
    val = [msg.values.flat[idx] for msg in msgs]
    grbs.close()
    return np.mean(val) if val else np.nan

def get_hrrr_series():
    """
    Download and extract 48h HRRR forecast at Boulder, CO.
    Returns DataFrame with time, wind speeds and directions at surface and 1000m, and surface temperature.
    """
    times, wind_sfc, wind_1000m, temp_sfc = [], [], [], []
    wind_dir_sfc, wind_dir_1000m = [], []
    for fh in range(1, HOURS+1):
        print(f"Processing HRRR fh={fh}")
        grib_path = download_hrrr_grib2(fh)
        # print(grib_path)
        if not grib_path:
            wind_sfc.append(np.nan)
            wind_1000m.append(np.nan)
            wind_dir_sfc.append(np.nan)
            wind_dir_1000m.append(np.nan)
            temp_sfc.append(np.nan)
            times.append(START + timedelta(hours=fh))
            continue
        # Surface winds and temp
        u_sfc = extract_hrrr_variable(grib_path, "10 metre U wind component", level=10)
        v_sfc = extract_hrrr_variable(grib_path, "10 metre V wind component", level=10)
        wind_sfc.append(np.hypot(u_sfc, v_sfc))
        # Wind direction (meteorological convention: direction FROM which the wind is blowing)  #!!! Need to add the direction difference between the HRRR grid and True North!
        dir_sfc = (180/np.pi) * np.arctan2(-u_sfc, -v_sfc)
        dir_sfc = (dir_sfc + 360) % 360
        wind_dir_sfc.append(dir_sfc)
        # 1000m winds (use 850hPa)
        u_1000 = extract_hrrr_variable(grib_path, "U component of wind", level=850)
        v_1000 = extract_hrrr_variable(grib_path, "V component of wind", level=850)
        wind_1000m.append(np.hypot(u_1000, v_1000))
        dir_1000 = (180/np.pi) * np.arctan2(-u_1000, -v_1000)
        dir_1000 = (dir_1000 + 360) % 360
        wind_dir_1000m.append(dir_1000)
        # Temp
        temp = extract_hrrr_variable(grib_path, "Temperature", level="surface")
        temp_sfc.append(temp-273.15)
        times.append(START + timedelta(hours=fh))
        os.remove(grib_path)
    return pd.DataFrame({
        "time": times,
        "hrrr_wind_sfc": wind_sfc,
        "hrrr_dir_sfc": wind_dir_sfc,
        "hrrr_wind_1000m": wind_1000m,
        "hrrr_dir_1000m": wind_dir_1000m,
        "hrrr_temp_sfc": temp_sfc
    })


### Google models


def init_ee():
    """
    Initialize Earth Engine.
    """
    if ee is None:
        raise ImportError("You must install earthengine-api: pip install earthengine-api")
    try:
        ee.Initialize(project=ee_project)   
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=ee_project)

def get_weathernext_series(model_type):
    """
    Download WeatherNExt (GraphCast or GeCat) series for Boulder, CO.
    model_type: GenCast or GraphCast
    var_dict: {output_name: gee_band}
    
        
    See documentation at :
        https://developers.google.com/earth-engine/datasets/catalog/projects_gcp-public-data-weathernext_assets_126478713_1_0#description
        https://developers.google.com/earth-engine/datasets/catalog/projects_gcp-public-data-weathernext_assets_59572747_4_0#description
    
    
    GraphCast: Forecast init times have 6 hour resolution (00z, 06z, 12z, 18z). Forecast lead times have 6 hour resolution up to a max lead time of 10 days.
    GenCast: Forecast init times have 6 hour resolution (00z, 06z, 12z, 18z). Forecast lead times have 12 hour resolution up to a max lead time of 15 days.
    
    graphcast_df = get_weathernext_series(collection_id = graphcast_id)
    gecat_df = get_weathernext_series(collection_id = gecat_id)
    
    """
    
    
    if model_type == "GraphCast":
        collection_id = "projects/gcp-public-data-weathernext/assets/59572747_4_0"  # GEE asset type
        forecast_hours = np.arange(0,HOURS,6) # GraphCast
    elif model_type == "GenCast":
        collection_id= "projects/gcp-public-data-weathernext/assets/126478713_1_0"  # GEE asset type
        forecast_hours = np.arange(0,HOURS,12) # GenCast
        #!!! Decide how to handle the ensemble members of GenCast. Currently just using one member. Make the average?
    else:
        print ("Model type not known, should be 'GraphCast' or 'GenCast'.")
        return

    var_dict = {  # Same variables/ bands for gencast and graphcast
        "u_10m": "10m_u_component_of_wind",  # 10m wind   
        "u_1000m": "1000_u_component_of_wind",  # wind at 1000m  
        "v_10m": "10m_v_component_of_wind",  # 10m wind   
        "v_1000m": "1000_v_component_of_wind",  # wind at 1000m   
        "temp_2m": "2m_temperature"  # 2m temp   
    }
     
    # Initalize google earth engine (GEE)
    init_ee()
    
    # Create a GEE geometry point at the coordinates (LON, LAT)
    point = ee.Geometry.Point([LON, LAT])  
    
    # For storing the results
    times = []
    results = {k: [] for k in var_dict}
    
    # Load an ImageCollection from Earth Engine specified by collection_id at START initialization time           
    if model_type == "GraphCast":
        imgcoll = ee.ImageCollection(collection_id)  \
                .filter(ee.Filter.date(START.isoformat()))\
                .filter(ee.Filter.bounds(point))
    elif model_type == "GenCast":
        imgcoll = ee.ImageCollection(collection_id)  \
                .filter(ee.Filter.date(START.isoformat()))\
                .filter(ee.Filter.bounds(point))\
                .filter(ee.Filter.eq('ensemble_member', '0'))   
                # .filter(ee.Filter.eq('forecast_hour', forecast_hours))  # !!! Works with one forecast hour, but not with multiple
                 
        
    # Convert the filtered collection into a list
    imglist = imgcoll.toList(imgcoll.size())   
    
    # Get the number of images in the collection 
    nimg = imgcoll.size().getInfo() 
    print (f"{nimg} images")
    
    
    # Loop over each image index i in the collection
    if nimg > 0:
        for i in range(nimg):
            print(i)
            try:
                # i-th image from imglist as an Earth Engine Image object
                img = ee.Image(imglist.get(i))
                
                # Image end_time info (end_time = valid time)     
                img_time = img.get('end_time').getInfo()
                
                                
                # # To get all property values of image
                # props_dict = img.toDictionary().getInfo()
                # print("All properties with values:")
                # for k, v in props_dict.items():
                #     print(f"{k}: {v}")
                
                # Temporary dict to store this image's band values
                vals_for_image = {}
                
                # Loop over each variable specified in var_dict as out_name and band (the band name in the image)
                for out_name, band in var_dict.items():
                    # Select the closest value to the point location - this is still necessary, because the image is still a raster, even after bounds filtering.
                    val = img.select(band).sample(
                        region=point,
                        scale=img.projection().nominalScale(),
                        numPixels=1,
                        geometries=True
                    ).first().getInfo()['properties']
                    
                    # Store result
                    vals_for_image[out_name] = val[band]
                
                # If all bands succeeded, append timestamp and band data
                times.append(img_time)
                for k, v in vals_for_image.items():
                    results[k].append(v)
                    
            except Exception as e:
                print(f"Skipping image {i} due to error: {e}")
            
            
    df = pd.DataFrame({"time": times, **results})
    
    # Calculate wind speed and direction
    df["wind_10m"] = np.hypot(df["u_10m"], df["v_10m"])
    df["wind_dir_10m"] =  (180/np.pi) * np.arctan2(-df["u_10m"], -df["v_10m"])
    df["wind_dir_10m"] = (df["wind_dir_10m"] + 360) % 360
    df["wind_1000m"] = np.hypot(df["u_1000m"], df["v_1000m"])
    df["wind_dir_1000m"] =  (180/np.pi) * np.arctan2(-df["u_1000m"], -df["v_1000m"])
    df["wind_dir_1000m"] = (df["wind_dir_1000m"] + 360) % 360 
    
    # Temperatrue in deg D
    df["temp_2m"] = df["temp_2m"] - 273.15
    
    df.time = pd.to_datetime(df.time).dt.tz_localize(None)
        
    return df 

# Tomorrow.io


def get_tomorrowio_series(lat, lon, start_time, hours, key):
    """
    Fetches a series of hourly temperature, wind speed, wind gust, and wind direction from Tomorrow.io API.

    Args:
        lat (float): Latitude
        lon (float): Longitude
        start_time (datetime): Start time (UTC)
        hours (int): Number of forecast hours

    Returns:
        pd.DataFrame: DataFrame with columns time, temp, wind_speed, wind_gust, wind_dir
    """
    from datetime import timedelta


    # Convert start_time to UTC ISO string
    start_iso = start_time.strftime("%Y-%m-%dT%H:%M:%SZ")
    end_time = start_time + timedelta(hours=hours)
    end_iso = end_time.strftime("%Y-%m-%dT%H:%M:%SZ")

    url = "https://api.tomorrow.io/v4/weather/forecast"
    params = {
        "location": f"{lat},{lon}",
        "fields": "temperature,windSpeed,windGust,windDirection",
        "timesteps": "1h",
        "startTime": start_iso,
        "endTime": end_iso,
        "apikey": TOMORROW_IO_API_KEY,
        "units": "metric"
    }
    
    # print URL: requests.Request('GET', url, params=params).prepare().url
    response = requests.get(url, params=params, verify=False)
    response.raise_for_status()
    data = response.json()

    # Parse the timeseries
    rows = []
    for interval in data.get("timelines", {}).get("hourly", []):
        values = interval.get("values", {})
        rows.append({
            "time": interval["time"],
            "tomorrowio_temp": values.get("temperature"),
            "tomorrowio_wind_speed": values.get("windSpeed"),
            "tomorrowio_wind_gust": values.get("windGust"),
            "tomorrowio_wind_dir": values.get("windDirection")
        })
      
    df = pd.DataFrame(rows)
    df['time'] = pd.to_datetime(df['time'])
    df['time'] = df['time'].dt.tz_localize(None)

    return df


### Meteomatics


def fetch_meteomatics(lat, lon, start, hours, username, password):
    """
    Fetch 10m wind, 1000m wind, and 2m temperature forecast from Meteomatics API.
    """
    # Time string for Meteomatics API: start--end:PT1H
    start_str = start.strftime("%Y-%m-%dT%H:%M:%SZ")
    end = start + timedelta(hours=hours-1)
    end_str = end.strftime("%Y-%m-%dT%H:%M:%SZ")
    time_str = f"{start_str}--{end_str}:PT1H"
    
    # Parameters
    params = [
        "wind_speed_10m:ms",          # 10m wind
        "wind_dir_10m:d",    # 10m wind dir
        "wind_gusts_10m_1h:ms",    # 10m wind gust 
        "t_2m:C"                      # 2m temperature
    ]
    param_str = ",".join(params)
    url = f"https://api.meteomatics.com/{time_str}/{param_str}/{lat},{lon}/json"
    
    response = requests.get(url, auth=(username, password), verify=False)
    if response.status_code != 200:
        print(f"Meteomatics API error: {response.status_code} {response.text}")
        return pd.DataFrame(columns=["time", "meteomatics_wind_10m", "meteomatics_wind_1000m", "meteomatics_temp_2m"])

    data = response.json()
    # Parse JSON to DataFrame
    times = [pd.to_datetime(t["date"]).tz_localize(None) for t in data["data"][0]["coordinates"][0]["dates"]]
    wind_10m = [t["value"] for t in data["data"][0]["coordinates"][0]["dates"]]
    wdir_10m = [t["value"] for t in data["data"][1]["coordinates"][0]["dates"]]
    wgust_10m = [t["value"] for t in data["data"][2]["coordinates"][0]["dates"]]
    temp_2m = [t["value"] for t in data["data"][3]["coordinates"][0]["dates"]]
    
    return pd.DataFrame({
        "time": times,
        "meteomatics_wind_10m": wind_10m,
        "meteomatics_wdir_10m": wdir_10m,
        "meteomatics_wgust_10m": wgust_10m,
        "meteomatics_temp_2m": temp_2m
    })





#%% Main program

# Fetch Meteomatics
print("Fetching Meteomatics 10m wind forecast...")
meteomatics_df = fetch_meteomatics(
    LAT, LON, START, HOURS, 
    username = meteomatics_user, 
    password = meteomatics_pw
)

# Fetch HRRR data
print("Fetching HRRR data...")
hrrr_df = get_hrrr_series()


# Fetch tomorrow.io
print("Fetching Tomorrow.io data...")
tomorrow_df = get_tomorrowio_series(lat = LAT, lon = LON, 
                           start_time = START, 
                           hours = HOURS,
                           key = TOMORROW_IO_API_KEY)

# Google models  

START = START - timedelta(hours = 6)   # !!! to be sure that data is available. Replace by last available data.

print("Fetching GraphCast (WeatherNext) data...")
graphcast_df = get_weathernext_series("GraphCast")  

print("Fetching GenCat (WeatherNext) data...")
gecat_df = get_weathernext_series(model_type = "GenCast")

START = START + timedelta(hours = 6)   # !!! remove later



### Merge all on time
df = hrrr_df.merge(meteomatics_df, on="time", how="outer")  \
            .merge(tomorrow_df, on="time", how="outer") \
            .merge(gecat_df.add_prefix("gencast_"), left_on="time", right_on = "gencast_time", how="outer") \
            .merge(graphcast_df.add_prefix("graphcast_"), left_on="time", right_on = "graphcast_time", how="outer")
df = df.sort_values("time").reset_index(drop=True)



### Get observation data
obs = read_synoptic(station_ids = [ "F7730"], synoptic_api_key = synoptic_api_key)  # "BLD02","F7730", "F1714", "ATC01","WXM8571" are other stations in Boulder



#%% Plot
fig, axs = plt.subplots(4, 1, figsize=(14, 9), sharex=True)


#!!! add initialization times seperately for each model

# --- Surface Wind ---
axs[0].plot(df["time"], df["hrrr_wind_sfc"], label="HRRR")
axs[0].plot(df["time"], df["meteomatics_wind_10m"], label="Meteomatics")
axs[0].plot(df["time"], df["tomorrowio_wind_speed"], label="Tomorrow.io")  
axs[0].plot(df["time"], df["graphcast_wind_10m"], ".", label="Graphcast")
axs[0].plot(df["time"], df["gencast_wind_10m"], ".", label="Gencast")
for station_id in obs.index.levels[0]:
    station_data = obs.loc[station_id]
    axs[0].plot(station_data.windspeed, "-", color = "black", label = "Obs " +station_data.iloc[0].stid)
axs[0].set_ylabel("Surface Wind (m/s)")
axs[0].legend()
axs[0].grid()
axs[0].axvline(START, color = "grey")

# --- 1000m Wind ---
axs[1].plot(df["time"], df["hrrr_wind_1000m"], label="HRRR")
axs[1].plot(np.nan, np.nan, label="")
axs[1].plot(np.nan, np.nan,  label="")  
axs[1].plot(df["time"], df["graphcast_wind_1000m"], ".", label="Graphcast")
axs[1].plot(df["time"], df["gencast_wind_1000m"], ".", label="Gencast")
axs[1].set_ylabel("1000m Wind (m/s)")
axs[1].legend()
axs[1].grid()
axs[1].axvline(START, color = "grey")

# --- Surface Temp ---
axs[2].plot(df["time"], df["hrrr_temp_sfc"], label="HRRR")
axs[2].plot(df["time"], df["meteomatics_temp_2m"], label="Meteomatics")
axs[2].plot(df["time"], df["tomorrowio_temp"], label="Tomorrow.io") 
axs[2].plot(df["time"], df["graphcast_temp_2m"], ".", label="Graphcast")
axs[2].plot(df["time"], df["gencast_temp_2m"], ".", label="Gencast")
for station_id in obs.index.levels[0]:
    station_data = obs.loc[station_id]
    axs[2].plot(station_data.temperature, "-", color = "black", label = "Obs " +station_data.iloc[0].stid)
axs[2].set_ylabel("Surface Temp (°C)")
axs[2].legend()
axs[2].grid()
axs[2].axvline(START, color = "grey")

# --- Wind Direction ---
axs[3].plot(df["time"], df["hrrr_dir_1000m"], ".", label="HRRR")
axs[3].plot(df["time"], df["meteomatics_wdir_10m"], ".", label="Meteomatics")
axs[3].plot(df["time"], df["tomorrowio_wind_dir"], ".", label="Tomorrow.io")  
axs[3].plot(df["time"], df["graphcast_wind_dir_10m"], ".", label="Graphcast")
axs[3].plot(df["time"], df["gencast_wind_dir_10m"], ".", label="Gencast")
for station_id in obs.index.levels[0]:
    station_data = obs.loc[station_id]
    axs[3].plot(station_data.winddirection, ".", color = "black", label = "Obs " +station_data.iloc[0].stid)
axs[3].set_ylabel("10m Wind Dir (°)")
axs[3].legend()
axs[3].grid()
axs[3].axvline(START, color = "grey")

axs[3].set_xlabel("Time (UTC)")
plt.suptitle("Boulder, CO: Model comparison")
plt.tight_layout()

