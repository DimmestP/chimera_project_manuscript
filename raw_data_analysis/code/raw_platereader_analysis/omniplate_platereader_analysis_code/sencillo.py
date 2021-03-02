
#!/usr/bin/env python
import sys, collections, re
try:
    import libsbml
    sbml= 1
except ImportError:
    sbml= 0


class sencillo:
    '''
    A python version of facile, which can be used at the command line with options:

    -M for a mathematica script for algebraic manipulation
    -y for python scripts to simulate ODEs with scipy
    -m for matlab scripts to simulate ODEs
    -j for julia scripts to simulate ODEs with DifferentialEquations
    -s for stochastic simulations with StochKit2
    -S for SBML version 2.3 and other versions as, for example, -S2.4
    -v show the version of Sencillo

    or within python with a typical work flow:

    from sencillo import sencillo
    s= sencillo('model.eqn')
    s.printmathematica()

    The attributes

    s.species (dictionary of molecular species with birth rate, death rate, and initial amount)
    s.alphspecies (alphabetical list of moecular species)
    s.variables (an ordered dictionary of variables and their values and modifiers)
    s.parameters (an ordered dictionary of parameters and their values)
    s.reactions (a dictionary of the individual reactions: type, name, products, rate constant, rate value, reactants)
    s.name (name of the .eqn file)

    are all created.

    If you use this software, please cite
        F Siso-Nadal, JF Ollivier, PS Swain. BMC Syst Biol 1:36 (2007)

    ---
    Specification of the model (format of the .eqn file)
    ---

    ---
    Chemical equations obeying mass action
    ---
    The .eqn has a similar format to the facile files (see the facile page at sourceforge http://facile.sourceforge.net/wiki/index.php/Main_Page). Examples of mass action reactions are

    binding: E + S <-> C ; 10; 2
    conversion: C -> P + E ; k= 1.2
    D -> null ; k
    null -> P ; 34

    ---
    Defining variables and parameters
    ---
    Parameters can be defined as constants that must not be functions of chemical species:

    parameter n= 300

    Variables can be functions of chemical species:

    variable gamma= gamma_max*E/(K_g + E)

    Note that variables can only be functions of species, parameters, and variables that have been defined previously.

    ---
    General chemical equations
    ---
    Mass action can be avoided by using => rather than -> and specifying reaction rates:

    C <=> F ; E*D*h ; j*n

    or

    variable r= u*F*E
    null => Z ; r
    A => B ; r

    if a variable is defined.

    ---
    Dimerization, etc.
    ---
    Reactions with multiple occurences of the same reactant should be written either explicitly or with an asterisk

    A + A -> A2 ; f
    C -> 4*A + E ; b
    n*A => B ; A/(A+K)

    ---
    Specifying initial conditions
    ---
    Initial conditions are specified by, for example,

    A= 1

    and any chemical species not specified will have an initial condition of zero.
    '''
    def __init__(self, fname= False):
        '''
        Parse file

        Arguments
        --
        fname: file to be opened
        '''
        self.ver= '1.03'
        if fname:
            self.name= fname.split('.')[0]
            self.reactions= {}
            self.species= {}
            parameters= []
            variables= []
            nononames= 0

            for line in open(fname):
                line= line.rstrip()
                if (not line or line == 'INIT' or line == 'EQN' or line[0] == '#'
                    or line[0] == '%'):
                    # ignore section headings and comments
                    continue
                else:
                    # remove any trailing semi-colons
                    if line[-1] == ';': line= line[:-1]
                    # remove any inline comments
                    if line.find('#') > -1: line= line[:line.find('#')]

                    # analyze chemical equation
                    if ';' in line:
                        # process chemical reaction
                        bits= line.split(';')

                        # assign reaction name
                        if ':' in bits[0]:
                            # reaction name is given
                            bitsofbits= bits[0].split(':')
                            reactionname= bitsofbits[0]
                            chemeq= bitsofbits[1].strip()
                        else:
                            # no reaction name: assign a numeric name
                            reactionname= str(nononames)
                            nononames += 1
                            chemeq= bits[0]

                        # find rates of reaction
                        chemrates= [crate.strip() for crate in bits[1:]]

                        # expand any asterisks, such as 2*A, in chemical reaction
                        if '*' in chemeq:
                            cmpts= []
                            for cb in chemeq.split():
                                if '*' in cb:
                                    # expand asterisk, for example to A + A
                                    cs= cb.split('*')
                                    cb= ' + '.join([cs[1],]*int(cs[0]))
                                cmpts.append(cb)
                            # make new expanded chemical equation
                            chemeq= ' '.join(cmpts)

                        # add reactions to dictionary of reactions
                        if '<->' in chemeq:
                            # reversible reaction: make two elementary reactions
                            chems= chemeq.split('<->')
                            chemeqf= chems[0] + '->' + chems[1]
                            chemeqb= chems[1] + '->' + chems[0]
                        elif '<=>' in chemeq:
                            # reversible reaction: make two elementary reactions
                            chems= chemeq.split('<=>')
                            chemeqf= chems[0] + '=>' + chems[1]
                            chemeqb= chems[1] + '=>' + chems[0]
                        if '<->' in chemeq or '<=>' in chemeq:
                            self.addreaction(reactionname + '_f', chemeqf, chemrates[0])
                            self.addreaction(reactionname + '_b', chemeqb, chemrates[1])
                        else:
                            # irreversible reaction
                            self.addreaction(reactionname, chemeq, chemrates[0])

                    else:
                        # not a reaction
                        bits= line.split()
                        if bits[0] == 'parameter':
                            # add parameter to list of parameters
                            parameters.append(' '.join(bits[1:]))
                        elif bits[0] == 'variable':
                            # add variables to list of variables
                            variables.append(' '.join(bits[1:]))
                        else:
                            # analyze initial conditions
                            bits= line.split('=')
                            self.species[bits[0].strip()]['initial amount']= bits[1].strip()

            # sort reactions so that rates with numerical values come first
            self.reactions= collections.OrderedDict(sorted(self.reactions.items(),
                                                      key= lambda x: x[1]['rate value']))

            # sort species so that initial values with numerical values come first
            self.species= collections.OrderedDict(sorted(self.species.items(),
                                                    key= lambda x: x[1]['initial amount']))

            # store alphabetically ordered species
            self.alphspecies= sorted(self.species.keys())

            # convert variables to a dictionary
            vardict= collections.OrderedDict()
            for var in variables:
                vardict[var.split('=')[0].strip()]= {'value' : var.split('=')[1].strip()}
            self.variables= vardict
            # add modifiers
            for var in vardict.keys():
                self.variables[var]['modifiers']= self.findmodifiers(vardict[var]['value'])

            # convert parameters to a dictionary
            pardict= collections.OrderedDict()
            for par in parameters:
                pardict[par.split('=')[0].strip()]= {'value' : par.split('=')[1].strip()}
            self.parameters= pardict



    #######
    # utility functions
    #######
    def addreaction(self, name, chemeq, chemrates):
        '''
        Takes a reaction and updates reactions and species dictionaries.

        Arguments
        --
        name: name of reaction, such as "degradation of A"
        chemeq: the chemical reaction, such as "A -> B"
        chemrates: the rate of the reaction, such as "k= 1"
        '''

        r, s= self.reactions, self.species
        if '->' in chemeq:
            bits= chemeq.split('->')
        else:
            bits= chemeq.split('=>')

        # find rate constant
        rbits= chemrates.split('=')
        if len(rbits) == 1:
            # no value given for rate constant
            rbits.append('')

        # find reactants and products
        reacs= [rsi.strip() for rsi in bits[0].split('+')]
        prods= [psi.strip() for psi in bits[1].split('+')]

       # add to dictionary of reactions
        r[name]= {'reactants' : reacs, 'products' : prods, 'rate value' : rbits[1].strip(),
                  'rate constant' : rbits[0].strip(), 'name' : chemeq}
        # birth reactions are special
        if 'null' in bits[0]: r[name]['reactants']= ''
        # death reactions are special
        if 'null' in bits[1]: r[name]['products']= ''

        # calculate mass action rate of reaction
        if '->' in chemeq:
            r[name]['mass action']= 1
            ratec= r[name]['rate constant']
            if '+' in ratec or '-' in ratec:
                # add brackets to avoid ambiguities
                ratec= '(' + ratec + ')'
            if 'null' in bits[0]:
                # birth reaction
                rate= ratec
            else:
                # first-order reaction
                rate= ratec + '*' + r[name]['reactants'][0]
                # higher-order reaction
                if len(r[name]['reactants']) > 1:
                    for species in r[name]['reactants'][1:]:
                        rate += '*' + species
        else:
            # the rate is specifed for => reactions
            r[name]['mass action']= 0
            if '+' in rbits[0] or '-' in rbits[0]:
                # add brackets to avoid ambiguities
                rate= '(' + rbits[0] + ')'
            else:
                rate= rbits[0]

        # add to dictionary of species
        for species in r[name]['reactants']:
            if species != '' and species != 'null':
                if species in s:
                    s[species]['death'] += '-' + rate
                else:
                    s[species]= {'birth' : '', 'death' : '-' + rate, 'initial amount' : '0'}
        for species in r[name]['products']:
            if species != '' and species != 'null':
                if species in s:
                    s[species]['birth'] += '+' + rate
                else:
                    s[species]= {'birth' : '+' + rate, 'death' : '', 'initial amount' : '0'}



    def findmodifiers(self, rate, exclude= []):
        '''
        Determines the species that are present in the expression of a chemical rate. These so-called modifiers are necessary for SBML. Note that some species may be indirectly acting as modifiers through variables.

        Arguments
        --
        rate: reaction rate
        exclude: these species are not modifiers
        '''
        bits= re.sub(r'[\^*\+-\/\(\)]', ';', rate).split(';')
        species= self.species.keys()
        varnames= self.variables.keys()
        mods= []
        for b in bits:
            bb= b.strip()
            if bb in species:
                mods.append(bb)
            elif bb in varnames:
                mods += self.variables[bb]['modifiers']
        # return unique list
        mods= list(set(mods))
        # carry out exclusions if required
        for rem in exclude: mods.remove(rem)
        return list(set(mods))



    def isfloat(self, s):
        try:
            float(s)
        except ValueError:
            return False
        return True




    #####
    def printmathematica(self, ofile= False):
        '''
        Output the model in Mathematica for analytical analysis (not simulation)

        Arguments
        --
        ofile: if specifed, output to this file
        '''
        s= self.species
        if ofile:
            out= open(ofile + '.ma', 'w')
        else:
            out= sys.stdout
        for sname in self.alphspecies:
            outstring= ('d' + sname + 'dt= ' + s[sname]['birth']
                        + s[sname]['death'] + ';\n')
            out.write(outstring)
        if ofile: out.close()


    #####
    def printodes(self, out, otype):
        '''
        Prints an ordinary differential equation version of the model to be called by a simulation routine.

        Arguments
        --
        out: destination for the output
        otype: either 'python' or 'matlab' or 'julia'
        '''
        p, r, s, v= self.parameters, self.reactions, self.species, self.variables
        if otype == 'python':
            deliml= '['
            delimr= ']'
            stp= ''
            inc= 0
        elif otype == 'matlab':
            deliml= '('
            delimr= ')'
            stp= ';'
            inc= 1
        elif otype== 'julia':
            deliml= '['
            delimr= ']'
            stp=''
            inc=1

        # define the vector of rates
        c= 0
        for rname in r.keys():
            # only define rates that have had values given
            if r[rname]['rate value'] != '':
                out.write('\t' + r[rname]['rate constant'] + '= rates' \
                          + deliml + str(c+inc) + delimr + stp +  '\n')
                c += 1
        out.write('\n')
                # define the vector of parameters
        for (c, par) in enumerate(p.keys()):
            out.write('\t' + par + '= parameters' + deliml \
                      + str(c+inc) +  delimr + stp + '\n')
        out.write('\n')
#            # define y as a vector of species
        for (c, sname) in enumerate(self.alphspecies):
            out.write('\t' + sname + '= y' + deliml + str(c+inc) + delimr + stp + '\n')
        out.write('\n')
            # print any variables
        for var in v.keys():
            if otype == 'python':
                out.write('\t' + var + '= ' + v[var]['value'].replace('^', '**') + stp + '\n')
            elif (otype == 'matlab') | (otype == 'julia'):
                out.write('\t' + var + '= ' + v[var]['value'] + stp + '\n')
                out.write('\n')
                # print the differential equations
        if otype == 'python':
            out.write('\tdydt= np.empty(' + str(len(s)) + ')\n')
        elif otype == 'matlab':
            out.write('\tdydt(size(y,1),1)= 0;\n')
        for (c, sname) in enumerate(self.alphspecies):
            if otype == 'python':
                out.write('\tdydt' + deliml + str(c+inc) + delimr + '= ' \
                          + s[sname]['birth'].replace('^', '**') \
                          + s[sname]['death'].replace('^', '**') + stp + '\n')
            elif (otype == 'matlab') | (otype == 'julia'):
                out.write('\tdydt' + deliml + str(c+inc) + delimr + '= ' \
                          + s[sname]['birth'] + s[sname]['death'] + stp + '\n')


    #####
    def printdriver(self, out, otype):
        '''
        Prints driver code for simulating an ordinary differential equation version of the model.

        Arguments
        --
        out: destination for the output
        otype: 'python', 'matlab' or 'julia'
        '''
        p, r, s, v= self.parameters, self.reactions, self.species, self.variables
        if otype == 'python':
            deliml= 'np.asarray(['
            delimr= '])'
            stp= ''
            fillc= ', '
            commt= '#'
        elif otype == 'matlab':
            deliml= '['
            delimr= ']'
            stp= ';'
            fillc= ' '
            commt= '%'
        if otype == 'julia':
            deliml= '['
            delimr= ']'
            stp= ''
            fillc= ', '
            commt= '#'

        # define parameters
        out.write(commt + ' parameters\n')
        outstring= ''
        for par in p.keys():
            out.write(par + '= ' + p[par]['value'] + stp + '\n')
            outstring += fillc + par
        out.write('parameters= ' + deliml + outstring[1:] + delimr + stp + '\n\n')

        # define rate constants
        out.write(commt + ' define rate constants\n')
        outstring= ''
        for rname in r.keys():
            # check that a rate value has been given before printing
            if r[rname]['rate value'] != '':
                out.write(r[rname]['rate constant'] + '= ' \
                          + r[rname]['rate value'] + stp + '\n')
                outstring += fillc + r[rname]['rate constant']
        out.write('rates= ' + deliml + outstring[1:] + delimr + stp + '\n\n')

        # define initial conditions
        out.write(commt + ' define initial conditions\n')
        outstring= '';
        for sname in self.alphspecies:
            if self.isfloat(s[sname]['initial amount']):
                out.write(sname + '_0= ' + s[sname]['initial amount'] + stp + '\n')
            else:
                if s[sname]['initial amount'] in self.alphspecies:
                    # uses a previously defined initial amount
                    out.write(sname + '_0= ' + s[sname]['initial amount'] + '_0' + stp + '\n')
                else:
                    out.write(sname + '_0= ' + s[sname]['initial amount'] + stp + '\n')
            outstring += fillc + sname + '_0'
        out.write('init= ' + deliml + outstring[1:] + delimr + stp + '\n\n')



    #####
    def printpython(self, ofile= False):
        '''
        Output the model in Python for simulation.
        A model file and a "driver" file that calls a simulation routine are created.

        Arguments
        --
        ofile: if specifed, output to this file
        '''
        s= self.species

        if ofile:
            out= open(ofile + '.py', 'w')
        else:
            out= sys.stdout
            ofile= self.name
        out.write('import numpy as np\n')
        out.write('import matplotlib.pylab as plt\n')
        out.write('from scipy.integrate import odeint\n\n')

        # function for odeint
        out.write('def ' + ofile + '(y, t, rates, parameters):\n')
        self.printodes(out, 'python')
        out.write('\n\treturn dydt\n')

        # driver code
        out.write('\n#######\n\n')
        self.printdriver(out, 'python')
        out.write('# call odeint (note args must be a tuple) \n')
        out.write('t= np.linspace(0, 100, 10)\n')
        out.write("species= ['" + self.alphspecies[0] + "'")
        for sname in self.alphspecies[1:]:
            out.write(", '" + sname + "'" )
        out.write(']\n')
        out.write('y= odeint(' + ofile + ', init, t, args= (rates, parameters), mxstep= 10000)\n')
        if not str(ofile): out.close()

    #####
    def printjulia(self, ofile= False):
        '''
        Output the model in Julia for simulation.
        A model file and a "driver" file that calls a simulation routine are created.

        Arguments
        --
        ofile: if specifed, output to this file
        '''
        s= self.species

        if ofile:
            out= open(ofile + '.jl', 'w')
        else:
            out= sys.stdout
            ofile= self.name
        out.write('using DifferentialEquations\n\n')

        # function for odeint
        out.write('function ' + ofile + '(dydt, y, (rates, parameters), t)\n')
        self.printodes(out, 'julia')
        # out.write('\n\treturn dydt\n')
        out.write('end\n')

        # driver code
        out.write('\n#######\n\n')
        self.printdriver(out, 'julia')
        out.write('# call solver (note args must be a tuple) \n')
        out.write('t= (0., 100.)\n')
        out.write("species= [:" + self.alphspecies[0])
        for sname in self.alphspecies[1:]:
            out.write(", :" + sname )
        out.write(']\n')

        out.write('prob= ODEProblem(' + ofile + ', init, t, (rates,parameters))\n')
        out.write('sol= solve(prob)')
        if not str(ofile): out.close()


    #####
    def printmatlab(self, ofile= False):
        '''
        Output the model in Matlab for simulation.
        A model file and a "driver" file that calls a simulation routine are created.

        Arguments
        --
        ofile: if specifed, output to this file
        '''
        s= self.species

        # function for ode23
        if ofile:
            out= open(ofile + '_odes.m', 'w')
        else:
            out= sys.stdout
            ofile= self.name
        out.write('function dydt= ' + ofile + '_odes(t, y, rates, parameters)\n\n')
        self.printodes(out, 'matlab')
        if str(ofile):
            print('\n\n')
        else:
            out.close()

        # driver file
        if not str(ofile): out= open(ofile + 'Driver.m', 'w')
        self.printdriver(out, 'matlab')
        out.write('% call solver routine \n')
        out.write('t0= 0;\n')
        out.write('tf= NO_TIME_SPECIFIED;\n')
        out.write('[t,y]= ode23s(@(t,y) ' + ofile + '_odes(t, y, rates, parameters), ' + '[t0 tf], init);\n')
        for (c, sname) in enumerate(self.alphspecies):
            out.write(sname + '= y(:,' + str(c+1) + ');\n')
        out.write('\n')
        if not str(ofile): out.close()




    #####
    def printstochkit(self, ofile= False):
        '''
        Output the model to StochKit for stochastic simulation.
        Variables cannot be used in the model description.

        Arguments
        --
        ofile: if specifed, output to this file
        '''
        #### Variables not implemented
        p, r, s, v= self.parameters, self.reactions, self.species, self.variables

        if ofile:
            out= open(ofile + '.xml', 'w')
        else:
            out= sys.stdout

        # add comment listing order of species
        out.write('<!-- ')
        porder= ' '
        for spec in sorted(s.keys()):
            porder += spec + ' '
        out.write(porder + ' -->\n')

        # preliminary details
        out.write('<Model>\n')
        out.write('\t<Description>' + self.name + '</Description>\n')
        out.write('\t<NumberOfReactions>' + str(len(r)) + '</NumberOfReactions>\n')
        out.write('\t<NumberOfSpecies>' + str(len(s)) + '</NumberOfSpecies>\n')

        # parameters
        out.write('\t<ParametersList>\n')
        if any(p):
            for par in p.keys():
                out.write('\t\t<Parameter>\n')
                out.write('\t\t <Id>' + par + '</Id>\n')
                out.write('\t\t <Expression>' + p[par]['value'] + '</Expression>\n')
                out.write('\t\t</Parameter>\n')
        for reac in r.keys():
            if r[reac]['rate value']:
                out.write('\t\t<Parameter>\n')
                out.write('\t\t <Id>' + r[reac]['rate constant'] + '</Id>\n')
                out.write('\t\t <Expression>' + r[reac]['rate value'] + '</Expression>\n')
                out.write('\t\t</Parameter>\n')
        out.write('\t</ParametersList>\n')

        # reactions
        out.write('\t<ReactionsList>\n')
        for reac in r.keys():
            out.write('\t\t<Reaction>\n')
            out.write('\t\t <Id>' + reac + '</Id>\n')
            out.write('\t\t <Description>' + r[reac]['name']+ '</Description>\n')
            # rate
            if r[reac]['mass action']:
                out.write('\t\t <Type>mass-action</Type>\n')
                out.write('\t\t <Rate>' + r[reac]['rate constant'] + '</Rate>\n')
            else:
                out.write('\t\t <Type>customized</Type>\n')
                out.write('\t\t <PropensityFunction>' + r[reac]['rate constant']
                          + '</PropensityFunction>\n')
            # reactants
            out.write('\t\t <Reactants>\n')
            if len(r[reac]['reactants']) == 1:
                    out.write('\t\t  <SpeciesReference id="' + r[reac]['reactants'][0] \
                              + '" stoichiometry="1"/>\n')
            elif len(r[reac]['reactants']) == 2:
                if r[reac]['reactants'][0] == r[reac]['reactants'][1]:
                    out.write('\t\t  <SpeciesReference id="' + r[reac]['reactants'][0]
                              + '" stoichiometry="2"/>\n')
                else:
                    out.write('\t\t  <SpeciesReference id="' + r[reac]['reactants'][0]
                          + '" stoichiometry="1"/>\n')
                    out.write('\t\t  <SpeciesReference id="' + r[reac]['reactants'][1]
                          + '" stoichiometry="1"/>\n')
            out.write('\t\t </Reactants>\n')
            # products
            out.write('\t\t <Products>\n')
            if len(r[reac]['products']) == 0:
                # degradation
                pass
            elif len(r[reac]['products']) == 1:
                out.write('\t\t  <SpeciesReference id="' + r[reac]['products'][0]
                    + '" stoichiometry="1"/>\n')
            else:
                if r[reac]['products'][0] == r[reac]['products'][1]:
                    if len(r[reac]['products']) > 2 and r[reac]['products'][0] == r[reac]['products'][2]:
                        out.write('\t\t  <SpeciesReference id="' + r[reac]['products'][0]
                                  + '" stoichiometry="3"/>\n')
                    else:
                        out.write('\t\t  <SpeciesReference id="' + r[reac]['products'][0]
                                  + '" stoichiometry="2"/>\n')
                else:
                    for i in range(len(r[reac]['products'])):
                        out.write('\t\t  <SpeciesReference id="' + r[reac]['products'][i]
                          + '" stoichiometry="1"/>\n')
            out.write('\t\t </Products>\n')
            out.write('\t\t</Reaction>\n')
        out.write('\t</ReactionsList>\n')

        # species
        out.write('\t<SpeciesList>\n')
        for spec in self.alphspecies:
            out.write('\t\t<Species>\n')
            out.write('\t\t <Id>' + spec + '</Id>\n')
            out.write('\t\t <InitialPopulation>' + s[spec]['initial amount']
                      + '</InitialPopulation>\n')
            out.write('\t\t</Species>\n')
        out.write('\t</SpeciesList>\n')

        out.write('</Model>\n')
        if ofile:
            out.close()

        # display order of species
        print('sencillo: Stochkit will display species alphabetically, i.e. as')
        print(porder)



    ####
    def printSBML(self, ofile, ver= ['2', '3']):
        '''
        Output the model to SBML.
        LibSBML must be installed.

        Arguments
        --
        ofile: if specifed, output to this file
        ver: version of SBML
        '''

        if sbml:

            p, r, s, v= self.parameters, self.reactions, self.species, self.variables

            # preliminaries: including SBML version, such as 2.3
            document= libsbml.SBMLDocument(eval(ver[0]), eval(ver[1]))
            model= document.createModel()
            c= model.createCompartment()
            c.setId('cell')
            c.setVolume(1.0)

            # create species
            for spec in self.alphspecies:
                ss= model.createSpecies()
                ss.setId(spec)
                ss.setCompartment('cell')
                ss.setInitialAmount(eval(s[spec]['initial amount']))

            # create parameters
            for par in p.keys():
                pp= model.createParameter()
                pp.setId(par)
                pp.setConstant(True)
                pp.setValue(eval(p[par]['value']))

            # variables must first be defined as parameters
            for var in v.keys():
                pp= model.createParameter()
                pp.setId(var)
                pp.setConstant(False)

            # define all rate constants as parameters
            for reac in r.keys():
                rvalue= r[reac]['rate value']
                if rvalue:
                    # a new rate constant is defined
                    newrate= r[reac]['rate constant']
                    pp= model.createParameter()
                    pp.setId(newrate)
                    if self.isfloat(rvalue):
                        # defined numerically
                        pp.setConstant(True)
                        pp.setValue(eval(rvalue))
                    else:
                        # defined algebraically
                        pp.setConstant(False)
                        vv= model.createAssignmentRule()
                        value= libsbml.parseFormula(rvalue)
                        vv.setVariable(newrate)
                        vv.setMath(value)

            # create variables
            for var in v.keys():
                vv= model.createAssignmentRule()
                value= libsbml.parseFormula(v[var]['value'])
                vv.setVariable(var)
                vv.setMath(value)

            # create reactions
            for reac in r.keys():
                rr= model.createReaction()
                rr.setId('r_' + reac)
                rr.setCompartment('cell')
                rr.setReversible(False)
                for i in range(len(r[reac]['reactants'])):
                    rrr= rr.createReactant()
                    rrr.setSpecies(r[reac]['reactants'][i])
                for i in range(len(r[reac]['products'])):
                    rrp= rr.createProduct()
                    rrp.setSpecies(r[reac]['products'][i])
                if r[reac]['mass action']:
                    # find rate by mass action
                    rrate= r[reac]['rate constant']
                    for reactant in r[reac]['reactants']:
                        rrate += '*'+reactant
                    mods= self.findmodifiers(rrate, exclude= list(set(r[reac]['reactants'])))
                else:
                    # rate already specified
                    rrate= r[reac]['rate constant']
                    mods= self.findmodifiers(rrate)
                # add modifiers
                for mod in mods:
                    mdf= rr.createModifier()
                    mdf.setSpecies(mod)
                rrk= rr.createKineticLaw()
                rrk.setFormula(rrate)

            # write file
            document.setModel(model)
            if ofile:
                libsbml.writeSBMLToFile(document, ofile + '.xml')
            else:
                libsbml.writeSBML(document, sys.stdout)
        else:

            print('sencillo: libSBML cannot be imported. SBML cannot be created.')






#####
if __name__ == '__main__':

    if len(sys.argv) == 3:

        fname= sys.argv[2]
        ofile= fname.split('.')[0]
        S= sencillo(fname)

        key= sys.argv[1].split('-')[1]
        if key == 'M':
            S.printmathematica(ofile)
        elif key == 'm':
            S.printmatlab(ofile)
        elif key ==  'y':
            S.printpython(ofile)
        elif key == 'j':
            S.printjulia(ofile)
        elif key ==  's':
            S.printstochkit(ofile)
        elif key[0] == 'S':
            if len(key) > 1:
                ver= key[1:].split('.')
            else:
                ver= ['2', '3']
            S.printSBML(ofile, ver)
        else:
            print(sencillo.__doc__)

    elif len(sys.argv) == 2 and sys.argv[1].split('-')[1] == 'v':
        print('Sencillo version number', sencillo().ver)

    else:

        print(sencillo.__doc__)
