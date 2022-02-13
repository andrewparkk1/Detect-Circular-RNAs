import re, os, stat

with open('sequence.txt', 'r') as rf:
        with open('100bp-circRNA.fa', 'w') as wf:
                while True:
                        #Organize line1
                        line = rf.readline()
                        if not line:
                                break;
                        arr = re.split(r'[-():>]', line)

                        #get the sequence
                        #negative strand
                        if arr[4] == "":
                                wf.write("> " + arr[1] + ":-:" + arr[2] + ":" + arr[3] + "\n")
                                line = rf.readline()
                                wf.write(line[:50].rstrip() + line[-51:] + "\n")
                        #positive strand
                        else:
                                wf.write("> " + arr[1] + ":" + arr[4] + ":" + arr[2] + ":" + arr[3] + "\n")
                                line = rf.readline()
                                wf.write(line[-51:].rstrip() + line[:50]+ "\n")


                        
