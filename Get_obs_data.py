import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import requests
from io import StringIO
import matplotlib.pyplot as plt



def read_metar(station, start_year, start_month, start_day, end_year, end_month, end_day):
    
    """
    use: 
        
    met = read_metar(station='TPH', 
                       start_year=year, start_month=month, start_day=day, 
                       end_year=day_after.year, end_month=day_after.month, end_day=day_after.day)
    
    year, month, date can be either strings or integers
    """
        
    import requests
    from io import StringIO
    
    # API endpoint (https://mesonet.agron.iastate.edu/request/download.phtml?network=NV_ASOS)
    url = "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py"
    params = {
        #'data': 'tmpf,dwpf,relh,drct,sknt,p01i,alti,mslp,vsby,gust,peak_wind_gust,peak_wind_drct,peak_wind_time',
        'data': 'drct,sknt,gust',
        'station': station,
        #'network': network,
        'latlon': "yes",
        'elev': "yes",
        'tz': 'UTC',
        # 'year1': start_year,
        # 'month1': start_month,
        # 'day1': start_day,
        # 'year2': end_year,
        # 'month2': end_month,
        # 'day2': end_day
    }
    
    # Send a GET request to the API
    response = requests.get(url, params=params)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Load the response content into a DataFrame
        csv_data = StringIO(response.text)
        metar = pd.read_csv(csv_data, delimiter=',', index_col=1, parse_dates=True, na_values="M")
    else:
        metar = pd.DataFrame()
        print(f"Failed to retrieve Metar data: {response.status_code}")
        
    # Ensure valid is a datetime index
    #metar = metar.drop_duplicates()

    # Wind speed in m/s
    metar["wspd"] =metar.sknt/1.94384    
    
    return metar




def read_synoptic(station_ids = ["BLD02","F7730"], synoptic_api_key = ""):
    
        
    #  https://api.synopticdata.com/v2/stations/timeseries?&token=efcf120a17ce45fba323a19de7b6a3f6&start=202508070000&end=202508080000&timeformat=%Y%m%d-%H%M&obtimezone=utc&units=metric&output=csv&stid=BLD02

    
    # Define the API URL
    url = "https://api.synopticdata.com/v2/stations/timeseries"
    
    # Initialize an empty list to store DataFrames
    combined_data = []

    
    for station in station_ids[:]: 
        
        print (station)
        
        
        params = {
            #   "country": "us",
            #  "state": "nm",
            "stid": station,
            # "complete": "1",
            "sensorvars": "1",
            #"start": "202404040000", # everything older than a year gets cut off.
            #"start": f"{year}01010000", # test with paid data access
            #"end":   "202504040100",
            #"end":   f"{year}12312359",
            "recent": "2880",  # in minutes
            "vars": "wind_speed,wind_direction,wind_gust,air_temp",
            "token": synoptic_api_key,
            "qc_checks " : "madis"
        }
    
        # Make the request
        response = requests.get(url, params=params, verify=False)
    
        # Check if successful
        if response.status_code == 200:
            
            try:
                data = response.json()
                
                # Save to file
                # with open(save_folder + f"/{station}_{year}.json", "w") as f:
                #     json.dump(data, f, indent=2)
                # print("Data saved")
                
                # Extract the station metadata and observation data
                station = data["STATION"][0]  # Assuming one station per file
                observations = station.get("OBSERVATIONS", {})
                
                # Extract metadata
                stid = station["STID"]
                name = station["NAME"]
                elevation = station["ELEVATION"]
                latitude = station["LATITUDE"]
                longitude = station["LONGITUDE"]
                state = station["STATE"]
                
                # Extract the position of the wind sensor (wind_speed_set_1)
                wind_speed_position = station["SENSOR_VARIABLES"]["wind_speed"]["wind_speed_set_1"].get("position")
                
                # Create a DataFrame from the observations
                df = pd.DataFrame(observations)
                
                # Remove '_set_1' from the column names
                df.columns = df.columns.str.replace('_set_1', '', regex=False)
                
                df.rename(columns={
                    'air_temp': 'temperature',
                    'wind_speed': 'windspeed',
                    'wind_direction': 'winddirection',
                    'wind_gust': 'gust'
                }, inplace=True)
                
                # Convert 'date_time' to datetime and set as index
                df['date_time'] = pd.to_datetime(df['date_time']).dt.tz_localize(None)
                df.set_index('date_time', inplace=True)
                
                # Add the station metadata as new columns
                df['stid'] = stid
                df['station_name'] = name
                df['elev'] = elevation
                df['lat'] = latitude
                df['lon'] = longitude
                df['state'] = state
                df['height'] = wind_speed_position
                
                df['elev'] = pd.to_numeric(df['elev'], errors='coerce') * 0.3048   # conversion from ft to m
                
                # Append the DataFrame to the list
                combined_data.append(df)
                
            except:
                continue
            
        else:
            print(f"Request failed with status code {response.status_code}")
            print(response.text)
            
            
            
    # Combine all DataFrames into a single DataFrame
    final_df = pd.concat(combined_data, ignore_index=False, keys = station_ids)
    

        
    # # Save the combined DataFrame to a CSV or Excel file
    # final_df.to_pickle('obs.pkl')
    
    return final_df


#%% Read and plot  data

if __name__ == "__main__":

    # Metar
    
    
    # metar = read_metar(station="BOCC2 ", 
    #                    start_year=2025, start_month=8, start_day=1, 
    #                    end_year=2025, end_month=8, end_day=2)
    # metar.wspd.plot()
    
    
    
    
    # Synoptic
    
    station_ids = ["BLD02","F7730", "F1714", "ATC01","WXM8571"]
    
    # check stations: https://viewer.synopticdata.com/map/data/now/air-temperature/F7730/about?layers=#map=11.44/40.0146/-105.2794
    
    obs = read_synoptic(station_ids)
    
    
    
    plt.figure()
    
    if not obs.empty and obs.index.nlevels > 1 and len(obs.index.levels[0]) > 0:
        for station_id in obs.index.levels[0]:
           
            station_data = obs.loc[station_id]
            
           # print(station_data)
            plt.plot(station_data.windspeed, "-", label = station_id)
    else:
        station_data = obs
        
    plt.legend()
        
    


    
    
    

