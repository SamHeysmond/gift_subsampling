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




