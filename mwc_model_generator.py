docstr = '''
  This code generates an MWC allosteric model of calmodulin
  in SBML format, using libSMBL. 
  
  Given the allosteric protein, it's states, binding sites, and
  allosteric targets, automatically generates a list of chemical
  species and association reactions.
  
  Output is a single .xml file.
  
  Massimo Lai - LeNovereLab, The Babraham Institute, Cambridge, UK

  Started on 2014-1-14
  
  2014/5/18 - Conformational transitions are now possible for target-
              bound molecules of calmodulin. (ML)
  
  -------------------------------------------
   I M P L E M E N T A T I O N    N O T E S:
  -------------------------------------------
  
  >> Each species' name is constructed according to the syntax:
     [Molecule]_[Conformation]_[Occupied_sites]_[Bound_target]

      - If no target is bound, the [Bound_target] field is '0'
      - If no ligand is bound, the [Occupied_sites] field is '0'
  
  >> Default ligand is calcium ('ca').  
  
  >> Species and reactions are generated automatically, as are
     kinetic parameter for CaM-target binding reactions.
  
  >> Manual editing is required for the  dictionaries containing "core 
     parameters" (with fixed value) and "secondary parameters" (i.e. 
     derived from the core parameters by means of assignment rules).
  
  >> IMPORTANT: as a modelling assumptions, all changes in affinity are
     assumed to be due to a change in koff. The "core parameters" for a 
     reactions are always Kd and kon, while koff is calculated as:
         
     koff = Kd * kon.
     
  >> Currently, information about parameter units is discarded 
     (no explicit units are given in the SBML file). Unit consistency 
     is the user's responsibility.
  
  
              
'''
import sys
import libsbml

XMLNAME = "wtcam_two_targets.xml" # edit to change output name

SBML_LEVEL = 2
SBML_VERSION = 4

ALLM       = 'cam' # name of the allosteric molecule
LIGAND     = 'ca'
NSITES     = 4
NLOBES     = 2
CONFSTATES = ['RR','RT','TR','TT']
LIGSTATES  = ['0','A','B','C','D','AB','AC','AD','BC',
              'BD','CD','ABC','ABD','ACD','BCD','ABCD']
# Edit this line to add targets; dummy target '0' must not be deleted
TARGETS    = ['0','ng','camkii'] 


# Indices (used to parse compound names) 
MINDEX = 0 # index for the "molecule" field 
CINDEX = 1 # index for the "conformation" field 
LINDEX = 2 # index for the "ligand" field 
TINDEX = 3 # index for the "target" field 

'''---------------------------------------------------------------------
  A L L O W E D   C O N F O R M A T I O N A L   T R A N S I T I O N S
  (we consider only transitions of individual lobes from T to R, and  
  assume that the reaction is reversible)

  TT == RT
  ||    ||
  TR == RR

  This dictionary describes the topology of the reactions that describe
  conformational transitions of the allosteric molecule.
  Each dict entry is in the form <current state> : <accessible states>
---------------------------------------------------------------------'''

CONF_TRANS = {
'TT' : ['RT','TR'],
'RT' : ['RR'],
'TR' : ['RR'],
}

'''---------------------------------------------------------------------
 C O R E   G L O B A L   P A R A M E T E R S   
 (from which all other parameters are calculated)
 each dictionary entry is in the format:    
 str(parameter) : ( float(value) , str(units) )
---------------------------------------------------------------------'''

CORE_PARAMETERS = {
'conc_cam' : (1e-6,None),
'conc_rbp' : (1e-6,None),
'conc_tbp' : (1e-6,None),
'kon_AR'   : (1e9,'per_mol_per_sec'),
'kon_AT'   : (1e9,'per_mol_per_sec'),
'kon_CR'   : (1e7,'per_mol_per_sec'),
'kon_CT'   : (1e7,'per_mol_per_sec'),
'kon_tbp'  : (1e8,'per_mol_per_sec'),
'kon_rbp'  : (1e8,'per_mol_per_sec'),
'KAT'      : (9.192e-5,'per_mol'),
'KBT'      : (9.192e-5,'per_mol'),
'KCT'      : (6.242e-5,'per_mol'),
'KDT'      : (6.242e-5,'per_mol'),
'Kd_tbp_RR': (1.0,'per_mol'),
'Kd_tbp_RT': (1.0,'per_mol'),
'Kd_tbp_TR': (1.0,'per_mol'),
'Kd_tbp_TT': (1.0,'per_mol'),
'Kd_rbp_RR': (1.0,'per_mol'),
'Kd_rbp_RT': (1.0,'per_mol'),
'Kd_rbp_TR': (1.0,'per_mol'),
'Kd_rbp_TT': (1.0,'per_mol'),
'lN'       : (398000,None),
'lC'       : (8616.61,None),
'cN'       : (0.000215,None),
'cC'       : (0.000317,None),
'k_R2T_N0'  : (10000,'per_sec'),
'k_R2T_C0'  : (10000,'per_sec'),
}
     
'''---------------------------------------------------------------------
 S E C O N D A R Y   P A R A M E T E R S    
 (formulas must be legal math expressions to be parsed by libSBML)    
 each dictionary entry is in the format:  str(parameter) : str(formula)
---------------------------------------------------------------------'''

SECONDARY_PARAMETERS = {
'KAR'        : 'KAT * cN',
'KBR'        : 'KBT * cN',
'KCR'        : 'KCT * cC',
'KDR'        : 'KDT * cC',
'kon_BR'     : 'kon_AR',
'kon_DR'     : 'kon_CR',
'kon_BT'     : 'kon_AT',
'kon_DT'     : 'kon_CT',
'koff_AR'    : 'KAR * kon_AR',
'koff_BR'    : 'KBR * kon_BR',
'koff_CR'    : 'KCR * kon_CR',
'koff_DR'    : 'KDR * kon_DR',
'koff_AT'    : 'KAT * kon_AT',
'koff_BT'    : 'KBT * kon_BT',
'koff_CT'    : 'KCT * kon_CT',
'koff_DT'    : 'KDT * kon_DT',
'koff_tbp_RR': 'Kd_tbp_RR * kon_tbp',
'koff_tbp_RT': 'Kd_tbp_RT * kon_tbp',
'koff_tbp_TR': 'Kd_tbp_TR * kon_tbp',
'koff_tbp_TT': 'Kd_tbp_TT * kon_tbp',
'koff_rbp_RR': 'Kd_rbp_RR * kon_rbp',
'koff_rbp_RT': 'Kd_rbp_RT * kon_rbp',
'koff_rbp_TR': 'Kd_rbp_TR * kon_rbp',
'koff_rbp_TT': 'Kd_rbp_TT * kon_rbp',
'k_R2T_N1'   : 'k_R2T_N0',
'k_R2T_C1'   : 'k_R2T_C0',
'k_R2T_N2'   : 'k_R2T_N0',
'k_R2T_C2'   : 'k_R2T_C0',
'k_T2R_N0'   : 'k_R2T_N0 / lN',
'k_T2R_C0'   : 'k_R2T_C0 / lC',
'k_T2R_N1'   : 'k_R2T_N1 / (lN * cN)',
'k_T2R_C1'   : 'k_R2T_C1 / (lC * cC)',
'k_T2R_N2'   : 'k_R2T_N2 / (lN * cN * cN)',
'k_T2R_C2'   : 'k_R2T_C2 / (lC * cC * cC)'
}
  

class CoreParameter:
    def __init__(self,name,value,units=None):
        self.name = name
        self.value = float(value)
        self.units = units
        
    def setName(self,name):
        self.name = name
        
    def setValue(self,value):
        self.value = value
        
    def setUnits(self,units):
        self.units = units
        
    def getName(self):
        return self.name
        
    def getValue(self):
        return self.value
        
    def getUnits(self):
        return self.units
        
    def hasUnits(self):
        if self.units == None: 
	    return False
        else: 
            return True 

class SecondaryParameter:
    def __init__(self,name,formula):
        self.name = name
        self.formula = formula
        
    def setName(self,name):
        self.name = name
        
    def getName(self):
        return self.name  
        
    def setFormula(self,formula):
        self.formula = formula
        
    def getFormula(self):
        return self.formula        
        
'''---------------------------------------------------------------------
 Class that stores the information  to describe a reversible binding 
 reaction between two chemical species, or monomolecular reaction with
 1 reactant. Mass action law is assumed, and  stoichiometry is 1 for all 
 reactants and products.
---------------------------------------------------------------------'''
class RevRxn:
    def __init__(self,reactants,products,kf,kb,label=''):
        self.reactants = list(reactants)
        self.products  = list(products)
        self.kf = str(kf)
        self.kb = str(kb)
        #print reactants,; sys.stdout.flush()
        #print products; sys.stdout.flush()
        reacstr = ' + '.join(self.reactants)
        prodstr = ' + '.join(self.products)
        self.name = reacstr + ' = ' + prodstr   
        if label == '':
            if len(reactants) == 2: 
                # Binding reactions            
                self.label='Binding of '+str(reactants[0])+' to '+ \
                str(reactants[1])
            elif len(reactants) == 1: 
                # Conformational transitions / Unimolecular reactions
                self.label='Transition from '+str(reactants[0])+' to '+\
                str(products[0])
        else:
            self.label = label  
            
    def __str__(self):
        s1 = self.label
        #s2 = ' + '.join(self.reactants) + ' = ' + ' + '.join(self.products)
        s2 = self.name        
        s3 = 'kon: %s,  koff: %s' %(self.kf, self.kb)        
        return '      '.join([s1,s2,s3])        
        
    def getName(self):
        return self.name           

    def buildRateLaw(self): # Mass Action Law
        fw = '(' + ' * '.join(self.reactants) + ')' + ' * ' + self.kf
        bw = '(' + ' * '.join(self.products)  + ')' + ' * ' + self.kb
        return fw + ' - ' + bw        
        
    def getReactants(self):
        return self.reactants

    def setReactants(self,reactants):
        self.reactants = reactants

    def getProducts(self):
        return self.products	

    def setProducts(self,products):
        self.products = products

    def getLabel(self):
        return self.label

    def setLabel(self,label):
        self.label = label

    
def addSpeciesToModel(model,comp,name,initAmount):
    
    sp = model.createSpecies()
    sp.setName(name)
    sp.setId(name)
    sp.setCompartment(comp.getId())
    sp.setInitialAmount(initAmount)
    
    return model

'''-------------------------------------------------------------------------
     Functions to help parse species names quickly
--------------------------------------------------------------------------'''

''' Return conformational state of the N lobe '''
def confN(species):
    fields = species.strip().split('_')
    conf   = fields[CINDEX]
    return conf[0]
    
''' Return conformational state of the N lobe '''   
def confC(species):
    fields = species.strip().split('_')
    conf   = fields[CINDEX]
    return conf[1]
    
''' Return conformational state of the N lobe '''   
def boundCalcium(species):
    fields = species.strip().split('_')
    ligand = fields[LINDEX]
    if ligand == '0':
        return 0.
    else: 
        return float(len(ligand))
        
''' Return conformational state of the N lobe '''          
def boundCalciumN(species):
    fields = species.strip().split('_')
    ligand = fields[LINDEX]
    nsites = ['A','B']
    ca = sum([s in nsites for s in ligand])
    return ca

def boundCalciumC(species):
    fields = species.strip().split('_')
    ligand = fields[LINDEX]
    csites = ['C','D']
    ca = sum([s in csites for s in ligand])
    return ca

''' Return the ligand CaM is bound to '''   
def boundTo(species,target):
    fields = species.strip().split('_')
    t = fields[TINDEX]
    if (target == t): return True
    else: return False

def unliganded(species):
    fields = species.strip().split('_')
    ligand = fields[LINDEX]
    target = fields[TINDEX]
    return ((ligand=='0') and (target=='0'))
    

    
'''-------------------------------------------------------------------------
     Functions to generate reactions and global quantities
--------------------------------------------------------------------------'''


''' Expand the global dictionaries of parameters to include affinities,
binding kinetics and conformational kinetics of the allosteric molecule
when bound to each target --- WARNING: QUICK AND DIRTY IMPLEMENTATION -- 
should be re-coded for generality'''
def expand_parameters():
    for t in TARGETS[1:]: # Skip the dummy target '0'
        # Generate core binding parameters
        kon = 'kon_'+t; CORE_PARAMETERS[kon] = (1e8,'per_mol_per_sec')
        Kd_RR = 'Kd_' + t + '_RR'; CORE_PARAMETERS[Kd_RR] = (1.0,'per_mol')
        Kd_RT = 'Kd_' + t + '_RT'; CORE_PARAMETERS[Kd_RT] = (1.0,'per_mol')
        Kd_TR = 'Kd_' + t + '_TR'; CORE_PARAMETERS[Kd_TR] = (1.0,'per_mol')
        Kd_TT = 'Kd_' + t + '_TT'; CORE_PARAMETERS[Kd_TT] = (1.0,'per_mol')
        # Generate secondary binding parameters:
        # 1 - dissociation constants
        koff_RR = 'koff_'+t+'_RR'; SECONDARY_PARAMETERS[koff_RR]=Kd_RR + '*' + kon
        koff_RT = 'koff_'+t+'_RT'; SECONDARY_PARAMETERS[koff_RT]=Kd_RT + '*' + kon
        koff_TR = 'koff_'+t+'_TR'; SECONDARY_PARAMETERS[koff_TR]=Kd_TR + '*' + kon
        koff_TT = 'koff_'+t+'_TT'; SECONDARY_PARAMETERS[koff_TT]=Kd_TT + '*' + kon
        # 1 - e-parameters for all transitions
        eNR = 'eNR_'+t;SECONDARY_PARAMETERS[eNR] = Kd_RR +'/'+ Kd_TR
        eNT = 'eNT_'+t;SECONDARY_PARAMETERS[eNT] = Kd_RT +'/'+ Kd_TT
        eCR = 'eCR_'+t;SECONDARY_PARAMETERS[eCR] = Kd_RR +'/'+ Kd_RT
        eCT = 'eCT_'+t;SECONDARY_PARAMETERS[eCT] = Kd_TR +'/'+ Kd_TT  
        # 2 - kinetic rate constants for conformational transitions, no ca2+ bound
        k_TT2TR_C0 = 'k_TT2TR_C0_'+t; SECONDARY_PARAMETERS[k_TT2TR_C0]='k_T2R_C0'+'/'+eCT
        k_RT2RR_C0 = 'k_RT2RR_C0_'+t; SECONDARY_PARAMETERS[k_RT2RR_C0]='k_T2R_C0'+'/'+eCR
        k_TT2RT_N0 = 'k_TT2RT_N0_'+t; SECONDARY_PARAMETERS[k_TT2RT_N0]='k_T2R_N0'+'/'+eNT
        k_TR2RR_N0 = 'k_TR2RR_N0_'+t; SECONDARY_PARAMETERS[k_TR2RR_N0]='k_T2R_N0'+'/'+eNR
        # 3 - kinetic rate constants for conformational transitions, 1 ca2+ bound
        k_TT2TR_C1 = 'k_TT2TR_C1_'+t; SECONDARY_PARAMETERS[k_TT2TR_C1]='k_T2R_C1'+'/'+eCT
        k_RT2RR_C1 = 'k_RT2RR_C1_'+t; SECONDARY_PARAMETERS[k_RT2RR_C1]='k_T2R_C1'+'/'+eCR
        k_TT2RT_N1 = 'k_TT2RT_N1_'+t; SECONDARY_PARAMETERS[k_TT2RT_N1]='k_T2R_N1'+'/'+eNT
        k_TR2RR_N1 = 'k_TR2RR_N1_'+t; SECONDARY_PARAMETERS[k_TR2RR_N1]='k_T2R_N1'+'/'+eNR
        # 4 - kinetic rate constants for conformational transitions, 2 ca2+ bound
        k_TT2TR_C2 = 'k_TT2TR_C2_'+t; SECONDARY_PARAMETERS[k_TT2TR_C2]='k_T2R_C2'+'/'+eCT
        k_RT2RR_C2 = 'k_RT2RR_C2_'+t; SECONDARY_PARAMETERS[k_RT2RR_C2]='k_T2R_C2'+'/'+eCR
        k_TT2RT_N2 = 'k_TT2RT_N2_'+t; SECONDARY_PARAMETERS[k_TT2RT_N2]='k_T2R_N2'+'/'+eNT
        k_TR2RR_N2 = 'k_TR2RR_N2_'+t; SECONDARY_PARAMETERS[k_TR2RR_N2]='k_T2R_N2'+'/'+eNR
        
    return   
        

''' Generate automatically the global quantities obtained by summation of a
large number of chemical species, and returns a list of strings that are the 
equations to be parsed:'''    
def generateGlobalQuantities(species):
    
    #create equation for total concentration of allosteric molecule:
    molname = ALLM    
    equations = []    
    
    # Generate formulas of total Ca sites, and Ybar of free CaM and CaM bound to each target,
    # and also for the single lobes    
    ydict  = {}
    yNdict = {}
    yCdict = {}

    RRdict = {}
    RTdict = {}
    TRdict = {}
    TTdict = {}

    totname = ALLM +'_tot'       
       
    # Generate formula for the total Ybar,YbarN and YbarC of all species:     
    key = 'ybar_tot'
    keyN ='ybarN_tot'
    keyC ='ybarC_tot'        
    boundCa = [boundCalcium(s) for s in species]
    boundCa_N =[boundCalciumN(s) for s in species]
    boundCa_C =[boundCalciumC(s) for s in species] 
    addends = [str(n)+'*'+l for l,n in zip(species,boundCa) if not n==0] #Total CaM
    addends_N = [str(n)+'*'+l for l,n in zip(species,boundCa_N) if not n==0] # N lobe only
    addends_C = [str(n)+'*'+l for l,n in zip(species,boundCa_C) if not n==0]   
    ydict[totname] = '(' + ' + '.join(species) + ')'
    yNdict[keyN] = '(' + ' + '.join(addends_N) + ') / ('+str(NSITES / NLOBES)+'*' + totname + ')'        
    yCdict[keyC] = '(' + ' + '.join(addends_C) + ') / ('+str(NSITES / NLOBES)+'*' + totname + ')'                
    ydict[key] = '(' + ' + '.join(addends) + ') / ('+str(NSITES)+'*' + totname + ')'
    
    # Generate formulas of RRbar, RTbar TRbar, TTbar for all species:
    keyRR = 'RRbar_' + 'tot'
    keyRT = 'RTbar_' + 'tot' 
    keyTR = 'TRbar_' + 'tot' 
    keyTT = 'TTbar_' + 'tot' 
    addends_RR = [s for s in species if confN(s)=='R' and confC(s)=='R']        
    addends_RT = [s for s in species if confN(s)=='R' and confC(s)=='T']
    addends_TR = [s for s in species if confN(s)=='T' and confC(s)=='R']
    addends_TT = [s for s in species if confN(s)=='T' and confC(s)=='T']    
    RRdict[keyRR] = '(' + ' + '.join(addends_RR) + ') / (' + totname +')'   
    RTdict[keyRT] = '(' + ' + '.join(addends_RT) + ') / (' + totname +')' 
    TRdict[keyTR] = '(' + ' + '.join(addends_TR) + ') / (' + totname +')' 
    TTdict[keyTT] = '(' + ' + '.join(addends_TT) + ') / (' + totname +')' 
    
    # For each target t:
    for t in TARGETS:
        listOfSpecies = [s for s in species if boundTo(s,t)]    
        # Generate formulas of Ybar, YbarN and YbarC of free CaM and CaM bound to target t:
        key = 'ybar_' + t
        keyN ='ybarN_'+ t
        keyC ='ybarC_'+ t
        boundCa = [boundCalcium(s) for s in listOfSpecies]
        boundCa_N =[boundCalciumN(s) for s in listOfSpecies]
        boundCa_C =[boundCalciumC(s) for s in listOfSpecies] 
        addends = [str(n)+'*'+l for l,n in zip(listOfSpecies,boundCa) if not n==0] #Total CaM
        addends_N = [str(n)+'*'+l for l,n in zip(listOfSpecies,boundCa_N) if not n==0] # N lobe only
        addends_C = [str(n)+'*'+l for l,n in zip(listOfSpecies,boundCa_C) if not n==0] # C lobe only
        btotname = ALLM + '_'+ t +'_tot'
        ydict[btotname] = '(' + ' + '.join(listOfSpecies) + ')'
        yNdict[keyN] = '(' + ' + '.join(addends_N) + ') / ('+str(NSITES / NLOBES)+'*' + btotname + ')'        
        yCdict[keyC] = '(' + ' + '.join(addends_C) + ') / ('+str(NSITES / NLOBES)+'*' + btotname + ')'                
        ydict[key] = '(' + ' + '.join(addends) + ') / ('+str(NSITES)+'*' + btotname + ')'
        
        # Generate formulas of RRbar, RTbar and TRbar, TTbar of free CaM and CaM bound to target t:
        keyRR = 'RRbar_' + t
        keyRT = 'RTbar_' + t 
        keyTR = 'TRbar_' + t 
        keyTT = 'TTbar_' + t 
        addends_RR = [s for s in listOfSpecies if confN(s)=='R' and confC(s)=='R']        
        addends_RT = [s for s in listOfSpecies if confN(s)=='R' and confC(s)=='T']
        addends_TR = [s for s in listOfSpecies if confN(s)=='T' and confC(s)=='R']
        addends_TT = [s for s in listOfSpecies if confN(s)=='T' and confC(s)=='T']   
        btotname = ALLM + '_'+ t +'_tot'
        RRdict[keyRR] = '(' + ' + '.join(addends_RR) + ') / (' + btotname +')'   
        RTdict[keyRT] = '(' + ' + '.join(addends_RT) + ') / (' + btotname +')' 
        TRdict[keyTR] = '(' + ' + '.join(addends_TR) + ') / (' + btotname +')' 
        TTdict[keyTT] = '(' + ' + '.join(addends_TT) + ') / (' + btotname +')' 
        
    # Generate formula for fraction of free CaM and CaM bound to each target:
    fdict = {}    
    for t in TARGETS:
        fkey = ALLM + '_' + t + '_bound_fraction'
        boundlist = [s for s in species if boundTo(s,t)]
        fdict[fkey] = '(' + ' + '.join(boundlist) + ') / ('+ totname + ')'
    
    # Now merge all the dictionaries containing addends for each formula:
    edict = dict(ydict.items() + yNdict.items() + yCdict.items() + \
                 RRdict.items() + RTdict.items()  + TRdict.items() + TTdict.items() + \
                 fdict.items())
    # Generate the equations list from the equation dictionary:
    equations += [k + ' = ' + edict[k] for k in edict.keys()]
    
    return equations
    
'''Generate all possible liganded and conformational configurations'''    
def generateSpecies(mol,confstates,ligstates,targets):
    return ['cam'+'_'+c+'_'+l+'_'+t for c in confstates \
             for l in ligstates for t in targets]

'''Generate reactions from each species, as the set of all 
possible reverible reactions that can produce each species'''
def generateBindingReactions(species):
    rxnlist = []
    for sp in species: 
        fields = sp.strip().split('_')
        mol    = fields[MINDEX]
        conf   = fields[CINDEX]
        ligand = fields[LINDEX]
        target = fields[TINDEX]
    	  # generate the reaction of association with the target
        if target != '0':
            newmol = '_'.join([mol,conf,ligand,'0'])
            # check that the generated species is legit:
            if newmol not in species:
                raise StandardError, 'generateReactions: error: illegal species'
            reactants = [target,newmol]
            products  = [sp]
            rxnlabel  = '%s binding to %s' %(target,newmol)
            # choose the correct kinetic parameters according to conformation:
            kf = 'kon_'+target
            kb = 'koff_'+target+'_'+conf
            rxn = RevRxn(reactants,products,kf,kb,rxnlabel)
            rxnlist.append(rxn)
        
        # generate the reaction of association with the ligand
        if ligand != '0':
            occsites = list(ligand)
            for site in ligand:
                newlig = occsites[:]
                newlig.remove(site)
                if len(newlig)==0: newlig = ['0']
                newmol = '_'.join([mol,conf,''.join(newlig),target])
                reactants = [LIGAND,newmol]
                products  = [sp]
                rxnlabel  = '%s binding to %s on site %s' %(LIGAND,newmol,site)
                # choose the kinetic param for  site and conformation:
                NCconf = list(conf)
                if site in ['A','B']: lobeconf = NCconf[0]
                elif site in ['C','D']: lobeconf = NCconf[1]
                kf = 'kon_'+site+lobeconf
                kb = 'koff_'+site+lobeconf
                rxn = RevRxn(reactants,products,kf,kb,rxnlabel)
                rxnlist.append(rxn)
        
    return rxnlist
    
#---------------------------------------------------------------------------    
''' Generate all reversible conformational transitions; to simplify 
    generation, we assume as direct reaction the transition T->R,
    and then make it reversible to include all the R->T transitions
    automatically.'''
def generateConformationalTransitions(species):
    ctrnlist = []    
    for s in species:
        if unliganded(s): 
            # For now let's only assume transitions in unliganded states
            fields = s.strip().split('_')
            mol    = fields[MINDEX]
            conf   = fields[CINDEX]
            ligand = fields[LINDEX]
            target = fields[TINDEX] 
            # Choose the new configuration based on accessible conformations:
            if conf in CONF_TRANS.keys():
                for newconf in CONF_TRANS[conf]:
                    reactant = [s]
                    product = ['_'.join([mol,newconf,ligand,target])]
                    # Check which lobe is changing conformation: 
                    if conf[0]==newconf[0]: 
                        lobe = 'C'             
                    elif conf[1]==newconf[1]: 
                        lobe = 'N'
                    else: 
                        print 'S O M E T H I NG \'S   W R O N G...'                    
                    # NB : simultaneous transition of 2 lobes is forbidden
                    kf = 'k_T2R_' + lobe + '0'
                    kb = 'k_R2T_' + lobe + '0'
                    ctrn = RevRxn(reactant,product,kf,kb)
                    ctrnlist.append(ctrn)
                   
        # Allow conformational transition in target-bound states
        if not unliganded(s): 
            # For now let's only assume transitions in unliganded states
            fields = s.strip().split('_')
            mol    = fields[MINDEX]
            conf   = fields[CINDEX]
            ligand = fields[LINDEX]
            target = fields[TINDEX]
            
            if target != '0': 
                # Choose the new configuration based on accessible conf:
                if conf in CONF_TRANS.keys():
                    for newconf in CONF_TRANS[conf]:
                        reactant = [s]
                        product = ['_'.join([mol,newconf,ligand,target])]
                        trname = conf+'2'+newconf
                        # Check which lobe is changing conformation: 
                        if conf[0]==newconf[0]: 
                            lobe = 'C'   
                            numca = boundCalciumC(s)
                        elif conf[1]==newconf[1]: 
                            lobe = 'N'
                            numca = boundCalciumN(s)
                        else: 
                            print 'S O M E T H I N G \'S   W R O N G...'                    
                        # NB: simultaneous transition of 2 lobes is forbidden
                        kf = 'k_' + trname + '_' + lobe + str(numca) + '_' + target
                        kb = 'k_R2T_' + lobe + str(numca)
                        ctrn = RevRxn(reactant,product,kf,kb)
                        ctrnlist.append(ctrn)
                
                
            else: #calcium but no target!
                # Choose the new configuration based on accessible conf:
                if conf in CONF_TRANS.keys():
                    for newconf in CONF_TRANS[conf]:
                        reactant = [s]
                        product = ['_'.join([mol,newconf,ligand,target])]
                        # Check which lobe is changing conformation: 
                        if conf[0]==newconf[0]: 
                            lobe = 'C'   
                            numca = boundCalciumC(s)
                        elif conf[1]==newconf[1]: 
                            lobe = 'N'
                            numca = boundCalciumN(s)
                        else: 
                            print 'S O M E T H I N G \'S   W R O N G...'                    
                        # NB: simultaneous transition of 2 lobes is forbidden
                        kf = 'k_T2R_' + lobe + str(numca)
                        kb = 'k_R2T_' + lobe + str(numca)
                        ctrn = RevRxn(reactant,product,kf,kb)
                        ctrnlist.append(ctrn)        
            
    return ctrnlist 
    
    
#---------------------------------------------------------------------------
''' Adds the core parameters to the model, stored in a dictionary
We call core parameters thos that cannot be derived by other parameters '''
def addCoreParameters(model,pardict):
    for p in pardict.keys():
        name = p
        value = float(pardict[p][0])
        units = str(pardict[p][1])
        asr = model.createParameter()
        asr.setId(name)
        asr.setValue(value)
#        if units != None: asr.setUnits(units)
    return model
    
    
#---------------------------------------------------------------------------    
''' Adds the secondary paramters, i.e. those parameters that are functions
of other parameters '''
def addSecondaryParameters(model,pardict):
    #First, create a list of global parameters and set their value to 0
    for p in pardict.keys():
        #First, create a global parameter and set its value to 0
        name = p
        formula = pardict[p]
        par = model.createParameter()
        par.setId(name)
        par.setValue(0)
        par.setConstant(False)
        #Then, create an assignment rule for each parameter:
        asr = model.createAssignmentRule()
        asr.setVariable(name)
        asr.setFormula(formula)
    return model
    

#---------------------------------------------------------------------------
def addGlobalQuantities(model,gqlist):
    for gq in gqlist:
        #create global parameter and set its value to 0:
        words = gq.split(' = ')
        quantity = words[0]
        equation = words[1]
        par = model.createParameter()
        par.setId(quantity)
        par.setValue(0)
        par.setConstant(False)
        #Then, create an assignment rule for each parameter:
        asr = model.createAssignmentRule()
        asr.setVariable(quantity)
        asr.setFormula(equation)
    return model
    
    
#-----------------------------------------------------------------------
def addReactionsToCompartment(model,rxnlist,comp):

    compName = comp.getId();

    for robj in rxnlist:           # Get the reaction data:
        name = str(robj.getName()); #print name
        rid = robj.getLabel(); #print rid
        reactants = robj.getReactants(); #print reactants
        products = robj.getProducts(); #print products
        kinlaw = robj.buildRateLaw(); #print kinlaw
        # Create the reaction in the model
        rxn = model.createReaction()
        rxn.setName(str(rid))
        rxn.setId('_'.join(rid.split())) #remove spaces
        for rct in reactants:
            react = rxn.createReactant()
            react.setSpecies(rct)
        for prd in products:
            prod = rxn.createProduct()
            prod.setSpecies(prd)
        # Multiply the given rate law by the compartment volume
        # this is needed for the SBML consistency check to succeed.
        formula = '(' + kinlaw + ') * ' + str(compName)
        # kinetic law
        kl = rxn.createKineticLaw()
        kl.setFormula(formula)
        
    return model
    
    
#---------------------------------------------------------------------------    
def addSpeciesToCompartment(species,model,compartment):
    compName = compartment.getId()
    for s in species:
        sp = model.createSpecies()
        sp.setId(s)
        sp.setName(s)
        sp.setCompartment(compName)
        sp.setInitialAmount(0)
    return model
    


#---------------------------------------------------------------------------    
outputfile='mwc_cam_hemiconcerted_rbp_tbp.xml'
sbmlDoc = libsbml.SBMLDocument(2, 4)
model = sbmlDoc.createModel()
model.setId('mwc_cam_hemiconcerted_rbp_tbp')

comp = model.createCompartment()
comp.setId('cytosol')
comp.setSize(1.0)   

#---------------------------------------------------------------------------
#        ---   D E F I N I T I O N   O F   U N I T S   ---
#---------------------------------------------------------------------------
# Not currently used
'''
#Create unit for koff parameters and conformational transitions:
unitdef =  model.createUnitDefinition()
unitdef.setId('per_second')
unit = unitdef.createUnit()
unit.setKind(libsbml.UNIT_KIND_SECOND)
unit.setExponent(-1)

# Create units for kon of bimolecular reactions:
unitdef =  model.createUnitDefinition()
unitdef.setId('litre_per_mole_per_second')
unit = unitdef.createUnit()
unit.setKind(libsbml.UNIT_KIND_LITRE)
unit.setExponent(1)
unit = unitdef.createUnit()
unit.setKind(libsbml.UNIT_KIND_MOLE)
unit.setExponent(-1)
unit = unitdef.createUnit()
unit.setKind(libsbml.UNIT_KIND_SECOND)
unit.setExponent(-1)

#Create unit for concentration and dissociation constants:
unitdef =  model.createUnitDefinition()
unitdef.setId('per_mol')
unit = unitdef.createUnit()
unit.setKind(libsbml.UNIT_KIND_MOLE)
unit.setExponent(1)
unit = unitdef.createUnit()
unit.setKind(libsbml.UNIT_KIND_LITRE)
unit.setExponent(-1)
'''
#---------------------------------------------------------------------------
#   ---  M A I N   ---
#---------------------------------------------------------------------------
expand_parameters()
camspecies = generateSpecies(ALLM,CONFSTATES,LIGSTATES,TARGETS)
rxns = generateBindingReactions(camspecies)
ctrs = generateConformationalTransitions(camspecies)
glqs = generateGlobalQuantities(camspecies)

allspecies = camspecies + [LIGAND] + TARGETS[1:]  # don't add the dummy '0' species
model = addSpeciesToCompartment(allspecies,model,comp)
model = addCoreParameters(model,CORE_PARAMETERS)
model = addSecondaryParameters(model,SECONDARY_PARAMETERS)
model = addReactionsToCompartment(model,rxns,comp) # reactions
model = addReactionsToCompartment(model,ctrs,comp) # conformational transitions
model = addGlobalQuantities(model,glqs)

#-----------------------------------------------------------------------
#   ---   C H E C K   M O D E L   A N D   W R I T E   S B M L   ---
#   
#  This code for consistency testing was taken from the libSBML website
#-----------------------------------------------------------------------
def validateExampleSBML(sbmlDoc):
    if (sbmlDoc == None):
        print("validateExampleSBML: given a None SBML Document");
        return False;

    consistencyMessages = "";
    validationMessages = "";
    noProblems = True;
    numCheckFailures = 0;
    numConsistencyErrors = 0;
    numConsistencyWarnings = 0;
    numValidationErrors = 0;
    numValidationWarnings = 0;

    # LibSBML 3.3 is lenient when generating models from scratch using the
    # API for creating objects.  Once the whole model is done and before it
    # gets written out, it's important to check that the whole model is in
    # fact complete, consistent and valid.

    numCheckFailures = sbmlDoc.checkInternalConsistency();
    if (numCheckFailures > 0):
        noProblems = False;
        for i in range(0, numCheckFailures):
            sbmlErr = sbmlDoc.getError(i);
            if (sbmlErr.isFatal() or sbmlErr.isError()):
                numConsistencyErrors = 1 + numConsistencyErrors;
            else:
                numConsistencyWarnings = 1+numConsistencyWarnings;
        sbmlDoc.printErrors();

    # If the internal checks fail, it makes little sense to attempt
    # further validation, because the model may be too compromised to
    # be properly interpreted.

    if (numConsistencyErrors > 0):
        consistencyMessages = consistencyMessages  + "Further validation aborted.";
    else:
        numCheckFailures = sbmlDoc.checkConsistency();
        if (numCheckFailures > 0):
            noProblems = False;
            for i in range(0, numCheckFailures):
                sbmlErr = sbmlDoc.getError(i);
                if (sbmlErr.isFatal() or sbmlErr.isError()):
                    numValidationErrors = 1+numValidationErrors;
                else:
                    numValidationWarnings = 1+numValidationWarnings;
            sbmlDoc.printErrors();

    if (noProblems):
        return True;
    else:
        if (numConsistencyErrors > 0):
			tmp = ""
			if (numConsistencyErrors > 1):
				tmp = "s";
			print("ERROR: encountered " + str(numConsistencyErrors) + " consistency error" + tmp + " in model '" + sbmlDoc.getModel().getId() + "'.");
		
        if (numConsistencyWarnings > 0):
			tmp = ""
			if (numConsistencyWarnings > 1):
				tmp = "s"
			print("Notice: encountered " + str(numConsistencyWarnings)
            + " consistency warning" + tmp
            + " in model '" + sbmlDoc.getModel().getId() + "'.");
        print(consistencyMessages);

        if (numValidationErrors > 0):
			tmp = ""
			if (numValidationErrors > 1):
				tmp = "s"
			print("ERROR: encountered " + str(numValidationErrors) + " validation error" + (tmp)
            + " in model '" + sbmlDoc.getModel().getId() + "'.");

        if (numValidationWarnings > 0):
			tmp = ""
			if (numValidationWarnings > 1):
				tmp = "s"
			print("Notice: encountered " + str(numValidationWarnings)
            + " validation warning" + (tmp)
            + " in model '" + sbmlDoc.getModel().getId() + "'.");
        print(validationMessages);

        return (numConsistencyErrors == 0 and numValidationErrors == 0);
# 
# 
# Writes the given SBMLDocument to the given file.
# 
# 
def writeExampleSBML(sbmlDoc, filename):
    result = libsbml.writeSBML(sbmlDoc, filename)
    if (result == 1):
        print("Wrote file '" + filename + "'")
        return True
    else:
        print("Failed to write '" + filename + "'")
        return False
        
       
SBMLok = validateExampleSBML(sbmlDoc);
if (SBMLok):
    writeExampleSBML(sbmlDoc, XMLNAME);
if (not SBMLok):
    print 'An error occurred -- SBML file not written'
