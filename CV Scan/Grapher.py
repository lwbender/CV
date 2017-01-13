import os
csvPath = os.path.expanduser('~') + "/Desktop/TEST1.csv"
dosId = '84C492E0DFCB'
timeArray = '09:22:39 AM,09:23:42 AM,09:24:19 AM,09:26:10 AM,'.split(',')
voltages = '0,1,2,3,'.split(',')



i = 0
import matplotlib.pyplot as plt
import datetime
import csv
inp = [] 
with open(csvPath) as file:
        for i in file:
                inp.append(i)
#print(len(inp))



timeArray = [datetime.datetime.strptime(datetime.datetime.today().strftime('%m/%d/%Y') + ' ' + x, '%m/%d/%Y %I:%M:%S %p') for x in timeArray[:-1]]
voltages = [float(y) for y in voltages[:-1]]

print('Times: ', timeArray)
print('Voltages: ', voltages)




PC1 = [] #i = 18
PC2 = [] #i = 19
PC3 = [] #i = 20
PC4 = [] #i = 21
PC5 = [] #i = 22
PC6 = [] #i = 23
PC7 = [] #i = 24
timeM = [] #i = -3


for i, q in enumerate(inp):
       # print(i, q)
        if i > 1:
                data = q.split(',')
                if data[6] == dosId:
                        PC1.append(float(data[18]))
                        PC2.append(float(data[19]))
                        PC3.append(float(data[20]))
                        PC4.append(float(data[21]))
                        PC5.append(float(data[22]))
                        PC6.append(float(data[23]))
                        PC7.append(float(data[24]))
                        timeM.append(datetime.datetime.fromtimestamp(float(data[-3])))



fullData = [PC1, PC2, PC3, PC4, PC5, PC6, PC7]
for q,i in enumerate(fullData):
        print("PC{0}".format(q))datetime.datetime.today()
        print(i)






pc1ByVolt = []
pc2ByVolt = []
pc3ByVolt = []
pc4ByVolt = []
pc5ByVolt = []
pc6ByVolt = []
pc7ByVolt = []



for ind, vTime in enumerate(timeArray):
        volt1 = []
        volt2 = []
        volt3 = []
        volt4 = []
        volt5 = []
        volt6 = []
        volt7 = []
        zeit = timeM[ind]

        print(vTime, zeit)
        
        if vTime < zeit:
                if ind + 1 < len(timeArray):
                        if zeit < timeArray[ind+1]:
                                print("Successfully appended data from index {0}".format(ind))
                                volt1.append(PC1[ind])
                                volt2.append(PC2[ind])
                                volt3.append(PC3[ind])
                                volt4.append(PC4[ind])
                                volt5.append(PC5[ind])
                                volt6.append(PC6[ind])
                                volt7.append(PC7[ind])
                else :
                        print("No further time")
                        volt1.append(PC1[ind])
                        volt2.append(PC2[ind])
                        volt3.append(PC3[ind])
                        volt4.append(PC4[ind])
                        volt5.append(PC5[ind])
                        volt6.append(PC6[ind])
                        volt7.append(PC7[ind]) 
        pc1ByVolt.append(volt1)
        pc2ByVolt.append(volt2)
        pc3ByVolt.append(volt3)
        pc4ByVolt.append(volt4)
        pc5ByVolt.append(volt5)
        pc6ByVolt.append(volt6)
        pc7ByVolt.append(volt7)

allDataByVolt = [pc1ByVolt, pc2ByVolt, pc3ByVolt, pc4ByVolt, pc5ByVolt, pc6ByVolt, pc7ByVolt]
#for debug in allDataByVolt:
#        print(debug)

allData = []
for pc in allDataByVolt:
        newPC = []
        for volt in pc:
                newPC.append(sum(volt)/len(volt))
        allData.append(newPC)





print(allData)
print(allData[0])








with open(os.path.expanduser('~') + "/Desktop/CVOutput" + datetime.datetime.now().strftime('%I%M%S') + '.csv', 'w', newline='') as csvOut:
        writer = csv.writer(csvOut)
        writer.writerow(["Measurement Time", "Voltage", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7"])
        for i, time in enumerate(timeM):
                print(i)
                writer.writerow([time, voltages[i], allData[0][i], allData[1][i], allData[2][i], allData[3][i], allData[4][i], allData[5][i], allData[6][i]])

for pc, i in enumerate(allData):
        plt.figure(pc+1)
        plt.plot(voltages, i)
        plt.xlabel("Voltage (V)")
        plt.ylabel("Capacitance (fF)")
        plt.title("PC{0}".format(pc+1))
        plt.grid(True)
plt.show(block=True)







	
