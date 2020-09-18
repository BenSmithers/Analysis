import os # OS operations
import pandas as pd # working with the excel file 

# path information to where the election results file is
file_dir = "/home/benito/Downloads"
file_name = "Election Test(1-3).xlsx"
path_to = os.path.join(file_dir, file_name)

if not os.path.exists(path_to):
    raise IOError("Could not find election results at {}".format(path_to))


data = pd.read_excel(path_to)

n_votes = len(data['ID'])

# retults will be stored as a dicitonary of dictionaries
# contains dict, keys for office
# each of those office dicts will have a dict with names as keys and values of votes 
results = {} 
logged_names = [] # stores a list of the names of voters
logged_email_prefix = [] # stores a list of email prefixes of voters 

skip_keys = ["ID", "Start time", "Completion time", "Email", "Name"]
offices = list(filter( lambda x : x not in skip_keys, data.keys()))

for vote in range(n_votes):
    flag = False

    # get the name
    if data["Name"][vote] in logged_names:
        flag = True
    else:
        logged_names.append(data["Name"][vote])

    # separate the domain from the email, log the username part
    email_sep = data["Email"][vote].split("@")[0]
    if email_sep in logged_email_prefix:
        flag = True
    else:
        logged_email_prefix.append( email_sep )

    # if the flag was set we print out a warning
    if flag:
        print("ATTENTION! It looks like someone voted twice")
        print("    Name: {}".format(data["Name"][vote]))
        print("   Email: {}".format(data["Email"][vote]))
        print("Their vote is not being counted")
        print()
        continue #skip to the next vote 

    # Log the choices for this voter's vote 
    for office in offices:
        if office not in results:
            results[office] = {}
        
        # Register each of the choices for office 
        choice = data[office][vote]
        if choice not in results[office]:
            results[office][choice] = 1
        else:
            results[office][choice] += 1
            
# print the results of the election
# very simple... 
for office in results:
    print("{}:".format(office))
    for choice in results[office]:
        print("    {} received {} votes".format(choice, results[office][choice]))

