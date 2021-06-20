import pandas as pd
import matplotlib.pyplot as plt
import crunchflow as cf


skip_lines = [0,1,2,3] # specify first few lines to be skipped


# vol = cf.tec('conc')
# print(vol.columns)

def col_title(file):
    # Instantiate empty list of columns
    columns = []

    # Open file and read line by line
    with open(file) as f:
        for line in f:
            if "TITLE" in line:
                title = line.split('"')[1]

            if "VARIABLES" in line:
                columns = line.split('"')

    # Remove x, y, z
    columns = [col.strip() for col in columns]
    for col in ['', 'VARIABLES =', 'X', 'Y', 'Z']:
        if col in columns:
            columns.remove(col)

    # Remove whitespace, which was reduced to len 0 via col.strip()
    remove_fields = [col for col in columns if len(col) == 0]

    for col in remove_fields:
        columns.remove(col)

    return columns, title

# C:\\Users\\nimaa\\OneDrive\\Desktop\\crunchflow\\workdesk\\
def parser(filename):
    with open('C:\\Users\\nimaa\\OneDrive\\Desktop\\crunchflow\\workdesk\\' + filename, 'r') as infile, open('result.csv', 'w') as outfile:
        for i,line in enumerate(infile):
            if i in skip_lines: continue
            outfile.write(" ".join(line.split()).replace(' ', ','))
            outfile.write("\n")
file_name = 'totcon4.dat'

col, title = col_title('C:\\Users\\nimaa\\OneDrive\\Desktop\\crunchflow\\workdesk\\' + file_name)
col = [c for c in col if c != ',']

parser(file_name)
cols = range(len(col)) # specify columns
x = pd.read_csv("result.csv", header = None, usecols=cols)

fig, ax = plt.subplots()
ax.plot(x.iloc[:,0], x.iloc[:,1:])
ax.set_xlabel('distance')
ax.set_ylabel('Concentration')
ax.legend(x)
plt.show()
plt.savefig('C:\\Users\\nimaa\\OneDrive\\Desktop\\crunchflow\\workdesk\\' + file_name[:-4], dpi = 80)