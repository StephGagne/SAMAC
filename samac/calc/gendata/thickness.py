################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####



#########################################################################        
#############################  thickness  ###############################
#########################################################################
def thickness(self):
    """ This property calculates the max and min cloud thickness (in meters) based on the height properties. The method defheight must be run at least once on the object for this property to work.
    The property can be accessed as CloudObj.thickness["min"] or CloudObj.thickness (returns a dict) """
    H=dict()
    try:
        H["max"]=list(); H["min"]=list();
        for i in range(len(self.props["height"])):
            H["max"].append(self.props["height"][i][3]-self.props["height"][i][0])
            H["min"].append(self.props["height"][i][2]-self.props["height"][i][1])
    except: print("[thickness] Height properties must be defined first using the defheight method.")
    return H
