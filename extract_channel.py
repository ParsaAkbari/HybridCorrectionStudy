import csv
import re
import commands
f = open('/home/pa354/c_elegans_data/analyse_nanocorr_data/accuracy_compare.csv', 'rb')
save_csv = open('/home/pa354/c_elegans_data/analyse_nanocorr_data/accuracy_channel.csv', 'w+')
csvwriter = csv.writer(save_csv, delimiter=',')
head = f.next().replace('\n', '').replace('"', '').split(" ")[1:]
head.append("channel_number")
csvwriter.writerow(head)
i = 1
for line in f:
    print(i)
    line = line.replace('\n', '').replace('"', "").split(" ")[1:]
    if "channel_" in line[1]:
        m = re.match("channel_(\d+).*", line[1])
        channel_number = m.group(1)+"_ebi"
    else:
        comm_output = commands.getstatusoutput("grep '^@"+line[1]+"' /home/pa354/minion_reads/tomas_longreads.fastq")[1]
        m = re.match(".*_ch(\d+)_.*", comm_output)
        channel_number = m.group(1)+"_miska"
    line.append(channel_number)
    csvwriter.writerow(line)
    i = 1 + i
