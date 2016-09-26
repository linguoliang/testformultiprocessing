with open("testdata.txt",'w') as outputfile:
    for x in range(20000000):
        outputfile.write(str(x)+'\n')