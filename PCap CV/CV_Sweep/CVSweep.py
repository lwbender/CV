import DDP_data as ddp
import veriipy
import veriipy.sql_tools as sqp
import veriipy.simple_filters as sfp
import veriipy.plotting_tools as ptp
import pytz
import csv
import matplotlib.pyplot as plt
import numpy as np
import time
from datetime import date, timedelta
import datetime


dosId = ''
with open('ini.txt') as f:
    for i in f:
        dosId = i
print("Accessing MySQL...")
Yesterday = date.today() - timedelta(days=1)
DDP = ddp.DDP_Data_sql(ble_id=dosId, start_date=Yesterday, end_date=None)
out = []
with open('o.csv') as f:
    reader = csv.reader(f,delimiter=',')
    for i in reader:
        out.append(i)
numM = len(out)
with open('CVSweep_data.csv', 'w', newline='') as newO:
    writer = csv.writer(newO)
    heading = ['MeasurementTimeStamp', 'ServerTimeStamp',
               'CVTimeStamp', 'Voltage',
               'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']
    writer.writerow(heading)
    sort_ddp = []
    for i, o in enumerate(out):
        newLine = []
        if i + 1 < len(out):
            tempDate0 = (str(date.today()) + " " + o[0]).strip()
            tempDate1 = (str(date.today()) + " " + out[i + 1][0]).strip()
            zeit = time.strptime(tempDate0, "%Y-%m-%d %I:%M:%S %p")
            zweit = time.strptime(tempDate1, "%Y-%m-%d %I:%M:%S %p")
            for w, ind in enumerate(DDP['timestamp']): 
                mal = time.strptime(str(ind), "%Y-%m-%dT%H:%M:%S.000000-0500")
                if mal > zeit and mal < zweit:
                    tempLine = [ind, DDP['creation_date'][w], o[0], float(o[1])/1000.0]
                    for x in range(1,7):
                        tempLine.append(DDP['PC' + str(x)][w])
                    newLine.append(tempLine)
                elif mal > zweit:
                    break
        for line in newLine:
            if line is not None or line != []:
                sort_ddp.append(line)
    for n in sort_ddp:
            writer.writerow(n)
   
    for x in range(1,7):
        if max([sort_ddp[y][3+x] for y in range(0, len(sort_ddp))]) > 100000:
            plt.figure(x)
            plt.plot([sort_ddp[z][3] for z in range(0, len(sort_ddp))], [sort_ddp[y][3+x] for y in range(0, len(sort_ddp))])
            plt.xlabel('Voltage (V)')
            plt.ylabel('Capacitance (fF)')
            plt.title('Capacitance v Voltage (PC' + str(x) + ')')
            plt.grid(True)
            plt.savefig("CVPC" + str(x) + "_Graph.png")
    plt.show(block=True)
