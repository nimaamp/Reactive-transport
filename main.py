import pandas as pd
import os
import plotly.express as px

skip_lines = [0,1,2,3]


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


    columns = [col.strip() for col in columns]
    for col in ['', 'VARIABLES =', 'X', 'Y', 'Z']:
        if col in columns:
            columns.remove(col)

    # Remove whitespace, which was reduced to len 0 via col.strip()
    remove_fields = [col for col in columns if len(col) == 0]

    for col in remove_fields:
        columns.remove(col)

    return columns, title


def parser(filename):
    with open('C:\\Users\\nimaa\\OneDrive\\Desktop\\crunchflow\\workdesk\\' + filename, 'r') as infile, open('result.csv', 'w') as outfile:
        for i,line in enumerate(infile):
            if i in skip_lines: continue
            outfile.write(" ".join(line.split()).replace(' ', ','))
            outfile.write("\n")




file_names = os.listdir('C:\\Users\\nimaa\\OneDrive\\Desktop\\crunchflow\\workdesk\\')
# Create Dictionary for File Name and Text
file_name_and_text = {}
for file in file_names:
    with open('C:\\Users\\nimaa\\OneDrive\\Desktop\\crunchflow\\workdesk\\' + file, "r") as target_file:
         file_name_and_text[file] = target_file.read()

    col, title = col_title('C:\\Users\\nimaa\\OneDrive\\Desktop\\crunchflow\\workdesk\\' + file)
    col = [c for c in col if c != ',']
    parser(file)
    cols = range(len(col))
    x = pd.read_csv("result.csv", header = None, usecols=cols)
    x.columns = ['Distance','H+','Cl-','Cs+', 'Na+' ,  'NO3-']


file_data = (pd.DataFrame.from_dict(file_name_and_text, orient='index')
             .reset_index().rename(index = str, columns = {'index': 'file_name', 0: 'text'}))

print(file_data)

x.columns = ['Distance','H+','Cl-','Cs+', 'Na+' ,  'NO3-']
#
# #Add new column for Animation process
# x['step'] = x.apply(lambda x: "one" , axis=1)
#
#
pd.options.plotting.backend = "plotly"
fig = px.line(x, x="Distance", y=['Distance','H+','Cl-','Cs+', 'Na+' ,  'NO3-'])


df1 = pd.DataFrame(x)
print(df1)
fig.show()










# file_name = 'totcon1.dat'
#
# col, title = col_title('C:\\Users\\nimaa\\OneDrive\\Desktop\\crunchflow\\workdesk\\' + file_name)
# col = [c for c in col if c != ',']
# parser(file_name)
# cols = range(len(col))
# x = pd.read_csv("result.csv", header = None, usecols=cols)
# x.columns = ['Distance','H+','Cl-','Cs+', 'Na+' ,  'NO3-']
#
# #Add new column for Animation process
# x['step'] = x.apply(lambda x: "one" , axis=1)
#
#
# pd.options.plotting.backend = "plotly"
# fig = px.line(x, x="Distance", y=['Distance','H+','Cl-','Cs+', 'Na+' ,  'NO3-'])
#
#
# df1 = pd.DataFrame(x)
# fig.show()
