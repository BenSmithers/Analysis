import pandas as pd

filename = r"/home/benito/Downloads/UTASalary.xlsx"
data = pd.read_excel(filename)

gta_names = ["GRADUATE TEACHING ASSISTANT I","GRADUATE TEACHING ASSISTANT II"]

n_entries = len(data) #5149 employees

def get_entries(data, **kwargs):
    """
    kwargs should be a (KEY IN DICT)=(VALUE OF KEY)
    
    returns a list of frames from data that match the presented keywords 
    """

    entries = []

    for i in range(n_entries):
        frame = pd.DataFrame(data, index=[i])

        for kwarg in kwargs: 
            if kwarg in frame:
                if not str(frame[kwarg][i]).lower() == kwargs[kwarg].lower():                    
                    continue
            else:
                continue
        entries.append(frame)
    return(entries)


ent = get_entries(data, Department="Physics")
print(ent[0])


