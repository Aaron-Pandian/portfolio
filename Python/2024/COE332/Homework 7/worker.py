#!/usr/bin/env python3

# Imports
from jobs import get_job_by_id, update_job_status, q, rd # Methods and clients
import redis
import time
from datetime import date, timedelta
import json

@q.worker
def do_work(jobid):
    update_job_status(jobid, 'in progress')
    
    # Initiate analysis
    job = get_job_by_id(jobid)
    start_date = job['start']
    end_date = job['end']
    incident_latitudes = []
    incident_longitudes = []

    # Gather the average location of all incidents witin timeframe
    for item in rd.keys():
        incident = json.loads(rd.get(item)) # Gets incident dictionary
        # Create check between timeframe
        data_date = incident['Published Date']
        incident_date = date(data_date[0:2], data_date[3:5], data_date[6:10])
        start_date = date(start_date[0:2], start_date[3:5], start_date[6:10]) # Month(2)/Day(2)/Year(4) add space at end
        end_date = date(end_date[0:2], end_date[3:5], end_date[6:10])
        if (start_date <= incident_date and incident_date <= end_date):
            # Append the locations to the lists
            incident_latitudes.append(float(incident['Latitude']))
            incident_longitudes.append(float(incident['Longitude']))
    
    summary = f"The average incident location is at ({sum(incident_latitudes/len(incident_latitudes))}N, {sum(incident_longitudes)/len(incident_longitudes)}W), and there were {len(incident_longitudes)} incidents during this time period."
    update_job_status(jobid, 'complete') 
    return summary

do_work()
