# THREADING PRACTISE
'''
import threading
import time

done =False

def worker(text):
    counter = 0
    while not done:
        time.sleep(1)
        counter+=1
        print(f"{text}: counter")

 # pass function as "object" not calling it -> we're passing the instance of it
 # daemon lets the script quit even if this thread is still running
 # args passed with a comma to be a TUPLE 
threading.Thread(target=worker, daemon=True, args=("ABC",)).start()
threading.Thread(target=worker, daemon=True, args=("XYZ",)).start()

input("press enter to quit")
done = True

##
# could also do.......

done =False

def worker(text):
    counter = 0
    while not done:
        time.sleep(1)
        counter+=1
        print(f"{text}: counter")

t1 = threading.Thread(target=worker, daemon=True, args=("ABC",)).start()
t2 =threading.Thread(target=worker, daemon=True, args=("XYZ",)).start()

t1.start()
t2.start()

t1.join()
t2.join()

# this will NOT reach this line of code because it only gets here once all the above threads have finished
input("press enter to quit")
done = True
'''

'''
## ASYNC PRACTISE CODE
# prints A,1,B,2
# NO CPU time is wasted
# take care when awaiting tasks
import asyncio

async def main():

    # this task is schedualed. when we have idle time it will call it
    task = asyncio.create_task(other_function())
    print("A")

    # since its sleeping here, it can go back and execute the task above
    await asyncio.sleep(1)

    print("B")

    # this needed to print 2 since other function goes to sleep and needs awaiting on?
    await task


async def other_function():
    print("1")
    await asyncio.sleep(2)
    print("2")

# need to call the function, not just refer to it so the brackets are needed
asyncio.run(main())
'''

# multiprocess practise
'''
import time
import multiprocessing

def sleep_fop_a_bit(seconds):
    print(f"Sleeping {seconds} second(s)")
    time.sleep(seconds)
    print("Done sleeping, am woke now")

# like cached variables? stores the info for the function/process
#p1 = multiprocessing.Process(target=sleep_fop_a_bit,args=[1])
#p2 = multiprocessing.Process(target=sleep_fop_a_bit,args=[1])

processes=[]

for x in range(10):
    p = multiprocessing.Process(target=sleep_fop_a_bit,args=[1])

    # multiprocessing requires this
    if __name__=='__main__':
        p.start()
        # pass in object
        processes.append(p)

for p in processes:
    p.join()

# exe time of the program up to this point
finish = time.perf_counter()
print("Finished running after seconds : ", finish)
'''

'''
# THREADING PRACTISE 2
import time
import concurrent.futures

arg1 ="SAM"
arg2="SMAS"
arg3="JOHN"

def speak(arg1,arg2,arg3):
    print(f"HELLO!! {arg1}")
    print(f"HELLO!! {arg2}")
    print(f"HELLO!! {arg3}")

img_names = [
    'photo-1516117172878-fd2c41f4a759.jpg',
    'photo-1532009324734-20a7a5813719.jpg',
    'photo-1524429656589-6633a470097c.jpg',
    'photo-1530224264768-7ff8c1789d79.jpg',
    'photo-1564135624576-c5c88640f235.jpg',
    'photo-1541698444083-023c97d3f4b6.jpg',
    'photo-1522364723953-452d3431c267.jpg',
    'photo-1513938709626-033611b8cc03.jpg',
    'photo-1507143550189-fed454f93097.jpg',
    'photo-1493976040374-85c8e12f0c0e.jpg',
    'photo-1504198453319-5ce911bafcde.jpg',
    'photo-1530122037265-a5f1f91d3b99.jpg',
    'photo-1516972810927-80185027ca84.jpg',
    'photo-1550439062-609e1531270e.jpg',
    'photo-1549692520-acc6669e2f0c.jpg'
]

t1 = time.perf_counter()


def process_image(img_name):

    img_name = "Broken_Image.png"
    speak(arg1,arg2,arg3)
    print(f'{img_name} was processed...')


with concurrent.futures.ThreadPoolExecutor() as executor:
    executor.map(process_image,img_names)


t2 = time.perf_counter()

print(f'Finished in {t2-t1} seconds')

'''
'''
import pandas
df = pandas.DataFrame()
df['Name'] = ['Anna', 'Pete', 'Tommy']
df['Scores'] = [97, 600, 200]
df['Questions'] = [2200, 75, 100]

print(df)

print("empty")
df= 0
del df
print(df)
'''

import concurrent.futures


PATH_TO_MAIN="/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/"
Total_GIFT_Mo98 = 0

Total_GIFT_Mo98_file=open(PATH_TO_MAIN+"Total_GIFT_Mo98.txt",'w')

Total_GIFT_Mo98_file.write(str(Total_GIFT_Mo98))

Total_GIFT_Mo98_file.close()

""" for _ in range(50):
    Total_GIFT_Mo98_file=open(PATH_TO_MAIN+"Total_GIFT_Mo98.txt",'r+')

    Total_GIFT_Mo98=Total_GIFT_Mo98_file.readlines()[-1]

    print("Total 1: ", Total_GIFT_Mo98)

    Total_GIFT_Mo98.replace("\n","")

    print("Total 2: ", Total_GIFT_Mo98)

    Total_GIFT_Mo98=int(Total_GIFT_Mo98)+1

    Total_GIFT_Mo98_file.write("\n")
    Total_GIFT_Mo98_file.write(str(Total_GIFT_Mo98))

    Total_GIFT_Mo98_file.close() """

""" numbers=[1,2,3,4,5,6,7,8,9,10]
shared_total=0

def process_num(number):

    Total_GIFT_Mo98_file=open(PATH_TO_MAIN+"Total_GIFT_Mo98.txt",'r+')

    Total_GIFT_Mo98=Total_GIFT_Mo98_file.readlines()[-1]

    shared_total+=1

    #print("Total 1: ", Total_GIFT_Mo98)

    Total_GIFT_Mo98.replace("\n","")

    #print("Total 2: ", Total_GIFT_Mo98)

    Total_GIFT_Mo98=int(Total_GIFT_Mo98)+1

    Total_GIFT_Mo98_file.write("\n")
    Total_GIFT_Mo98_file.write(str(Total_GIFT_Mo98))

    Total_GIFT_Mo98_file.close()

    print(f"number {number} done with!")



with concurrent.futures.ThreadPoolExecutor() as executor:
    executor.map(process_num,numbers,shared_total)

print("finished threading")
 """

# THIS WORKED!

""" #import modin.pandas as pd
import numpy as np
import pandas
import sys
import modin

# Implement your answer here. You are also free to play with the size
# and shape of the DataFrame, but beware of exceeding your memory!

import modin.pandas as pd

frame_data = np.random.randint(0, 100, size=(2**10, 2**5))
df = pd.DataFrame(frame_data)

# ***** Do not change the code below! It verifies that 
# ***** the exercise has been done correctly. *****

try:
    assert df is not None
    assert frame_data is not None
    assert isinstance(frame_data, np.ndarray)
except:
    raise AssertionError("Don't change too much of the original code!")
assert "modin.pandas" in sys.modules, "Not quite correct. Remember the single line of code change (See above)"

import modin.pandas
assert pd == modin.pandas, "Remember the single line of code change (See above)"
assert hasattr(df, "_query_compiler"), "Make sure that `df` is a modin.pandas DataFrame."

print("Success! You only need to change one line of code!") """

""" 
# works?
import modin.pandas as pd
import pandas

#############################################
### For the purpose of timing comparisons ###
#############################################
import time
import ray
# Look at the Ray documentation with respect to the Ray configuration suited to you most.
ray.init()
#############################################


start = time.time()

pandas_df = pandas.read_csv("yellow_tripdata_2015-01.csv", parse_dates=["tpep_pickup_datetime", "tpep_dropoff_datetime"], quoting=3)

end = time.time()
pandas_duration = end - start
print("Time to read with pandas: {} seconds".format(round(pandas_duration, 3)))

##########################################
start = time.time()

modin_df = pd.read_csv("yellow_tripdata_2015-01.csv", parse_dates=["tpep_pickup_datetime", "tpep_dropoff_datetime"], quoting=3)

end = time.time()
modin_duration = end - start
print("Time to read with Modin: {} seconds".format(round(modin_duration, 3)))

print("Modin is {}x faster than pandas at `read_csv`!".format(round(pandas_duration / modin_duration, 2)))



start = time.time()

big_pandas_df = pandas.concat([pandas_df for _ in range(25)])

end = time.time()
pandas_duration = end - start
print("Time to concat with pandas: {} seconds".format(round(pandas_duration, 3)))

start = time.time()

big_modin_df = pd.concat([modin_df for _ in range(25)])

end = time.time()
modin_duration = end - start
print("Time to concat with Modin: {} seconds".format(round(modin_duration, 3)))

print("Modin is {}x faster than pandas at `concat`!".format(round(pandas_duration / modin_duration, 2))) """

import modin.pandas as pandas
import ray


ray.init(_plasma_directory="/tmp") # setting to disable out of core in Ray
""" mypath="/mnt/c/users/sheys/OneDrive/LIFE_4137_ProjectDS/Project_DS_Files/Github/gift_subsampling/"
this_df=pandas.read_csv(mypath+"test_text.csv",sep=",",na_filter=False,usecols=["Age"])
print(this_df.head())
this_df.to_csv(mypath+"written_test_text.csv",header=True,index=False)
 """
print("before")

this_df=pandas.DataFrame(columns=["A","B"])
this_df.insert(1,"REPLACED B","REPLACED")
print(this_df.head())
this_df["SUBSAMPLE_NUM"] = 999

print("after")
print(this_df.head())
print("length")
print(len(this_df))