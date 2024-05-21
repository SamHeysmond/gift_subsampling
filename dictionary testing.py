# dictionary testing
dictionary_1= {
    101230: 0.001,
    123313: 0.002,
    123341: 0.00002
}

test_file=open("dictionary_test_output","w")

for position,p_value in dictionary_1.items():
    print(str(position)+","+str(p_value)+"\n")
    test_file.write(str(position)+","+str(p_value)+"\n")

# test sorting it out
sorted_dict={}
for key in sorted(dictionary_1,key=dictionary_1.get):
    sorted_dict[key]=dictionary_1[key]

print("Before sort: ", dictionary_1)
print("AFTER sort: ", sorted_dict)

dictionary_1=sorted_dict
print("AFTER sort: ", dictionary_1)

print("end of script")