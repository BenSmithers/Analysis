import numpy as np
import sys
import pandas as pd 

class Vote:
    """
    This class represents a person's vote for a category. It is constructed with this ordered list 

    It contains a ranked list of their votes, and a list of candidates which are no longer in the running. 
    When the 'get' function is called, it returns the highest-ranked choice that hasn't been dropped. 
    When the 'drop' function is called, it adds a new candidate to the list of candidates no longer running 
    """
    def __init__(self, votes):
        if not isinstance(votes, (list, tuple)):
            raise TypeError("Expected votes to be {}, not {}".format(list, type(votes)))

        # copy the ordered list of votes. Python's weird with references, so we do this to be careful 
        self.ordered_votes = [str(vote) for vote in votes]
        self.skip = []

    def drop(self, key):
        if key not in self.skip:
            self.skip.append(key)

    def get_all(self):
        """
        Returns a list of candidates which haven't been dropped 
        """
        votes = []
        for each in self.ordered_votes:
            if each not in self.skip:
                votes.append(each)
        return(votes)

    def get(self):
        index = 0
        while self.ordered_votes[index] in self.skip:
            index += 1 

            if index==len(self.ordered_votes):
#                print(self.ordered_votes)
#                print(self.skip)
                return

        return(self.ordered_votes[index])
    
    def __str__(self):
        return("Rank-{} Vote for {}".format(len(self.ordered_votes), self.get()))

    def __repr__(self):
        return(self.__str__())

"""
Columns have names like "Best Music [First Chice]"
So we have these two utility functions to separate the rank of the category and the name of the category. 
We do keep the right bracket (]), but that's okay. 
"""

def extract_category( column_title ):
    """
    Gets the name of the category from an entry in the header
    """
    return( column_title.split("[")[0])

def extract_rank( column_title ):
    """
    Gets the rank of this column (first, second, etc...)
    """
    return( column_title.split("[")[1])

def get_winner(ballot):
    """
    The ballot entry should be a dictionary with entries for each candidate, and values being the number of votes. 

    It gets the sum of all the votes, and checks if any candidate has a **majority** of votes. Returns the key if so, returns None if there is no winner 
    """
    if not isinstance(ballot, dict):
        raise TypeError("Expected {}, not {}".format(dict, type(ballot)))

    total = sum(ballot.values())
    for entry in ballot.keys():
        if float(ballot[entry])/total > 0.5:
            return entry


def build_category_dict(header, dmode=False):
    """
    I use this 
    """
    cats = {}
    for entry in header:
        if entry=='Timestamp' or entry=="Email Address":
            continue
        
        new_cat = extract_category(entry) # entry.split("[")[0]
        cats[new_cat] = {} if dmode else []
    return(cats)

def assign_votes( submission, header, cats):
    """
    This function parses a row from the excel file. It builds and assigns a person's votes for each category 

    submission - the person's list of votes (raw data!)
    header - a list of categories
    cats - the datastructure containing all the votes 

    Note! This assumes the submission's choices are ordered First->Second->Third->Fourth (etc.) 
    Will produce inaccurate results if this is not the case.
            Improvement: extract the rank, parse it to a number, and insert it in the appropriate place on the list 
    """

    last_key = "" # we look for a change in the category name to know when a Vote object is ready to be made and submitted 
    wip_vote = []
    for column in range(len(submission)):
        if column==0:
            continue
        if column==1:
            continue
       
        
        key = extract_category(header[column]) 
        if key!=last_key: # new category, meaning we should wrap up, build the vote, and add it to the cats object 
            if last_key!="":
                build_vote = Vote(wip_vote)
                cats[last_key].append(build_vote)
            
                wip_vote = []
        else:
            pass

        # Add this column to the WIP vote. 
        wip_vote.append(submission[column])
        last_key = key
    
    # need to make sure that last entry is added! 
    build_vote = Vote(wip_vote)
    cats[key].append(build_vote)
    
    return(cats)    

def extract_vote_new(categories, vote_obj):
    last_key = ""
    wip_vote = []
    print(vote_obj.keys())

    for column in vote_obj.keys():
        if column=="Timestamp" or column=="Email Address":
            continue
        key = extract_category(column)
        if key!=last_key and last_key!="":
            categories[last_key].append(Vote(wip_vote))
            wip_vote = []
        wip_vote.append(vote_obj[column])
        last_key = key
    categories[last_key].append(Vote(wip_vote))

    return categories



def load(filename="tb2023.xlsx"):
    """
    This function loads the file and prints out the winners 
    """

   # data = np.loadtxt(filename, dtype=str, delimiter=",")
    data = pd.read_excel(filename)
    
    ##header = data[0]
    header = data.keys()
    

    # first we parse the data, building a dictionary with keys for each category - and whose values are lists of Vote objects 
    categories = build_category_dict(header)
    data = data.transpose()

    if False:
        skipped_header = False
        for row in data:
            if not skipped_header:
                skipped_header = True
                continue

            categories = assign_votes( row, header, categories) 
    
    for index in data:
        categories = extract_vote_new(categories, data[index])

    at_least_one_bad = False

    print("Ranked Choice Voting")
    for category in categories.keys():
        bad = False # flag for "something bad happened"

        votes = categories[category]
        if len(votes)==0:
            bad = True
            print("Category {} had no votes".format(category))
            continue

        rnd = 0     # loop counter 
        while True:
            ballot = {} # used to keep track of the number of votes for each candidate 
            
            # we need to make sure we have an entry (even if zero) for each candidate still in the running 
            for candidate in votes[0].get_all():
                ballot[candidate] = 0

            # assign the votes
            for vote in votes: 
                choice = vote.get()
                if choice is None:
                    # this happens when a vote no longer has a valid option
                    # occurs whenever there is a tie. Doesn't need a warning
                    bad = True
                    continue
                else:
                    if choice in ballot:
                        ballot[choice] += 1
                    else:
                        ballot[choice] = 1
            if bad:
                break 

            # check for a winner. 
            winner = get_winner(ballot)
            if winner is not None: 
                break
            else:
                rnd +=1 
                if rnd>=7:
                    print("Category {} exceeded loop counter!".format(category))
                    bad = True
                    break
            
            # if we're here, there was no winner. So now we drop the worst-performing candidate(s)
            min_val = min(ballot.values())
            for key in ballot.keys():
                if ballot[key]<=min_val:
                    for vote in votes:
                        vote.drop(key)
        
        # print the results
        if not bad:
            print("{}:    {}".format(category, get_winner(ballot)))
        else:
            print("Couldn't find winner for {}".format(category))
            at_least_one_bad = True
   
    # this is the backup vote-couting system for when ranked choice fails
    if at_least_one_bad:
        print()
        print("Starting SIMPLE Vote System")
        skipped_header = False
        bcat = build_category_dict(header, True)
        for row in data:
            if not skipped_header:
                skipped_header = True
                continue
            for col in range(len(row)):
                if col==0:
                    continue # skip timestamp
                movie = row[col]
                cat = extract_category(header[col])
                raw_value = extract_rank(header[col]) 
                value = 0
                
                mapping = ["first", "second","third","fourth","fifth","sixth", "seventh", "eighth", "ninth", "tenth"]
                for entry in range(len(mapping)):
                    if mapping[entry] in raw_value.lower():
                        value = -1*(1+entry)
                        break
                    if entry == (len(mapping)-1):
                        print("Didn't find a mapping...")
                

                if movie not in bcat[cat]:
                    bcat[cat][movie] = value
                else:
                    bcat[cat][movie] += value 
        for cat in bcat.keys():
            winner = None
            winning_val = None
            print("Category: {}".format(cat))
            print(bcat[cat])


load()
