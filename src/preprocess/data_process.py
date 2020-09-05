import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import operator
import matplotlib.lines as mlines
import pylab as pl

class data_process:

    def __init__(self):
        self.df_students = []
        self.df_nominal = []
        self.activity_mapping = {'Advance catheter':'z', 'Anatomic identification':'l', 'Anesthetize': 'o', 
                        'Blood return': 'q','Check catheter position': '2', 'Check flow and reflow': '1',
                        'Check wire in long axis': 'w', 'Check wire in short axis':'v',
                        'Clean puncture area':'d', 'Compression identification':'n', 'Cover probe':'h',
                        'Doppler identification': 'm', 'Drap puncture area':'e',
                        'Drop probe': 'r', 'Gel in probe': 'g','Get in sterile clothes': 'c',
                        'Guidewire install':'t', 'Hand washing':'b','Position patient': 'k',
                        'Position probe':'j', 'Prepare implements':'a', 'Puncture':'p',
                        'Put sterile gel':'i', 'Remove guidewire':'0', 'Remove syringe':'s',
                        'Remove trocar':'u', 'Ultrasound configuration': 'f',
                        'Widen pathway': 'y', 'Wire in good position':'x'}
        self.map_activity_stage = {}
    # load CSV
    def load_data(self, student_file='CCC19 - Log CSV.csv',nominal_file='CCC19_Log_normative_simulated_v6.csv'):
        self.df_students = pd.read_csv(student_file, error_bad_lines=False,delimiter=',')
        self.df_nominal = pd.read_csv(nominal_file, error_bad_lines=False,delimiter=',')
    
    def df_students_process(self):
        if self.df_students.empty != True:
            df_students = self.df_students
            
            df_students['ACTIVITY_short'] = df_students.ACTIVITY.map(self.activity_mapping)
            df_activity_stage = df_students[['ACTIVITY','ACTIVITY_short', 'STAGE']]
            df_activity_stage = df_activity_stage.groupby(['ACTIVITY','ACTIVITY_short','STAGE']).size().reset_index().rename(columns={0:'count'})
            df_activity_stage.sort_values(['STAGE'])
            student_activities = df_students['ACTIVITY'].unique()
            self.stage_list = df_students['STAGE'].unique()
            # create a dict to map activites to Stages
            
            keys = df_activity_stage.ACTIVITY.values
            values = df_activity_stage.STAGE.values
            for idx, key in enumerate(keys,start=0):
                self.map_activity_stage[key] = values[idx]
            #print(df_students)	
            df_students_by_case = df_students.groupby('CASEID').aggregate(lambda x: list(x))
            # add and edit some columns 
            numOfEvents = [] # per case
            start_datetime = [] 
            end_datetime = []
            time_diff = [] # time taken for each event
            time_total = [] # total time per case
            Ressource = [] 
            Round = []
            for idx, date_list in enumerate(df_students_by_case.START, start = 0):
                Ressource.append(df_students_by_case.RESOURCE[idx][0])
                Round.append(df_students_by_case.ROUND[idx][0])
                start_datetime.append([datetime.datetime.strptime(x, '%m/%d/%Y %H:%M:%S') for x in date_list])
                end_datetime.append([datetime.datetime.strptime(x, '%m/%d/%Y %H:%M:%S') for x in df_students_by_case.END[idx]])
                time_diff.append([end_datetime[idx][t] - start_datetime[idx][t] 
                                for t, k in enumerate(df_students_by_case.START[idx], start = 0)])
                time_total.append(sum(time_diff[idx][:],datetime.timedelta()))
                numOfEvents.append(len(df_students_by_case.ACTIVITY[idx]))
            df_students_by_case.START = start_datetime
            df_students_by_case.END = end_datetime
            df_students_by_case['delta_t'] = time_diff
            df_students_by_case['total_time'] = time_total
            df_students_by_case['numOfEvents'] = numOfEvents
            df_students_by_case.RESOURCE = Ressource
            df_students_by_case.ROUND = Round
            student_naming = {'R_13_1C': '1', 'R_14_1D': '2', 'R_21_1F': '3',
                'R_31_1G': '4', 'R_32_1H': '5', 'R_33_1L': '6',
                'R_45_2A': '7', 'R_46_2B': '8', 'R_47_2C': '9',
                'R_48_2D': '10'}
            df_students_by_case.RESOURCE = df_students_by_case.RESOURCE.map(student_naming)
            #create empty lists for each stage
    
            stage_lists_dict = dict()
            for stage in self.stage_list:
                stage_lists_dict[stage] = {'activities':[],'time': []}
    
            # for each case create a list of activities for each stage --> list of lists --> new column
            VIDEOEND = [t for t in df_students_by_case.VIDEOEND]
            VIDEOSTART =[t for t in df_students_by_case.VIDEOSTART]
            VIDEOTIME = [list(np.array(VIDEOEND[idx]) - np.array(VIDEOSTART[idx])) for idx in range(len(VIDEOSTART))]
            df_students_by_case['VIDEOTIME'] = VIDEOTIME    
            for idx, row in df_students_by_case.iterrows():
                for stage in self.stage_list:
                    index_stage = [index for index, value in enumerate(row.STAGE) if value == stage]
                    stage_lists_dict[stage]['activities'].append(row.ACTIVITY_short[min(index_stage):max(index_stage)+1])
                    stage_lists_dict[stage]['time'].append(row.VIDEOTIME[min(index_stage):max(index_stage)+1])
                    
            for stage in self.stage_list:
                df_students_by_case[stage] = stage_lists_dict[stage]['activities']
                df_students_by_case[stage+'_time'] = stage_lists_dict[stage]['time']
            
            # make new index/ case ID
            df_students_by_case['CaseID'] = df_students_by_case['RESOURCE'].str.cat(df_students_by_case['ROUND'],sep="_")
            df_students_by_case['CaseID']
            df_students_by_case.set_index(keys = ['CaseID'], inplace = True)
            df_students_by_case.RESOURCE = df_students_by_case.RESOURCE.astype(int)
            df_students_by_case = df_students_by_case.sort_values(['RESOURCE','ROUND'])  
            
            return 	df_students_by_case
            
        else:
            print("Please load the student's CSV file correctly")
            
        
    def df_nominal_process(self):
        if self.df_nominal.empty != True:
            df_nominal = self.df_nominal
            nominal_activities = df_nominal['event'].unique()
            # delete invisible activities
            df_nominal = df_nominal[(df_nominal['event'] != 'INVISIBLE No Return') & (df_nominal['event'] != 'INVISIBLE No good position')]
            nominal_activities = df_nominal['event'].unique()
            student_activities = self.df_students['ACTIVITY'].unique()
            np.sort(student_activities) ==  np.sort(nominal_activities)
            df_nominal['stage'] = df_nominal['event']
            df_nominal.stage = df_nominal.stage.map(self.map_activity_stage)
    
            df_nominal['ACTIVITY_short'] = df_nominal.event.map(self.activity_mapping)
            df_normal_filtered = df_nominal.drop(columns = ['simulated:logProbability','concept:simulated'])
            # groupby
            df_normal_filtered_by_case = df_normal_filtered.groupby('case')
            # activity lists
            df_normal_filtered_by_case = df_normal_filtered.groupby('case').aggregate(lambda x: list(x))
            # get number of events per case
            numOfEvents = []
            for idx, event_list in enumerate(df_normal_filtered_by_case.event, start = 0):
                numOfEvents.append(len(event_list))
            df_normal_filtered_by_case['numOfEvents'] = numOfEvents
            #create empty lists for each stage
    
            case_complete_list = []
            # for each case create a list of activities for each stage --> list of lists --> new column
                
            for idx, row in df_normal_filtered_by_case.iterrows():
                if row.event[-1] == 'Check catheter position':
                    case_complete_list.append(idx)
            #print('completed cases in sumulation: ', len(case_complete_list))
            df_normal_filtered_by_case_completed = df_normal_filtered_by_case.iloc[case_complete_list,:].reset_index()
            #create empty lists for each stage
    
            stage_lists_dict = dict()
            for stage in self.stage_list:
                stage_lists_dict[stage] = []
    
            # for each case create a list of activities for each stage --> list of lists --> new column
                
            for idx, row in df_normal_filtered_by_case_completed.iterrows():
                for stage in self.stage_list:
                    index_stage = [index for index, value in enumerate(row.stage) if value == stage]
                    try:
                        stage_lists_dict[stage].append(row.ACTIVITY_short[min(index_stage):max(index_stage)+1])
                    except:
                        print('failed for row: ',idx, ' and stage: ', stage )			
            # add new columns to df
            for stage in self.stage_list:
                df_normal_filtered_by_case_completed[stage] = stage_lists_dict[stage]
            
            ##remove duplication
            activity_str = []
            for l in df_normal_filtered_by_case_completed.ACTIVITY_short:
                str1 = ''.join(str(e) for e in l)
                activity_str.append(str1)
            activity_set = list(set(activity_str))
            activity_index = []
            for i in activity_set:
                activity_index.append(activity_str.index(i))
            
            activity_index.sort()
            df_normal_filtered_by_case_completed=df_normal_filtered_by_case_completed.iloc[activity_index].reset_index()
            df_normal_filtered_by_case_completed = df_normal_filtered_by_case_completed.drop('index', 1)
            
            return df_normal_filtered_by_case_completed
            
        else:
            print("Please load the nominal CSV file correctly")	
            
            
			