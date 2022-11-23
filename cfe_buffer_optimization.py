from opentrons import protocol_api
from doepy import build
from math import ceil
import pandas as pd

metadata = {
    'apiLevel': '2.8',
    'protocolName': 'Cell-free expression Buffer optimization setup',
    'description': '''This protocol is for setting up a range 
    of different buffer compositions to optimize the buffer
    solution for the cell-free expression. It takes the premix and adds
    of Mg-glutamate, K-glutamate and PEG8000. In total 343 combinations are tested.''',
    'author': 'Karen Therkelsen (s173684@dtu.dk)',
}

nsamples = 343 + 1

def run(protocol):
    
    ## Load instrument, modules and labware ##

    #--------------#--------------#--------------#
    # 2mLEppendorf #   384-well   #    Trash     #
    #  tubes, 4C   #   plate      #              #
    #--------------#--------------#--------------#
    #  PCR strips, #   P20 tips   #  P300 tips   #
    #  4C          #              #              #
    #------------- #------------- #------------- #
    #              # 15mL falcon  #              #
    #              # rack         #              #
    #------------- #------------- #------------- #
    #              #              #              #
    #              #              #              #
    #------------- #------------- #------------- #


    # Pipettes and tips
    tips20 = protocol.load_labware('opentrons_96_tiprack_20ul', 8)
    tips300 = protocol.load_labware('opentrons_96_tiprack_300ul', 9)
    p20 = protocol.load_instrument('p20_single_gen2', mount='left', tip_racks=[tips20])
    p300 = protocol.load_instrument('p300_single_gen2', mount='right', tip_racks=[tips300])


    # 384-well plate
    plate = protocol.load_labware('corning_384_wellplate_112ul_flat', 11)
    
    # Tube rack
    rack = protocol.load_labware ('opentrons_15_tuberack_falcon_15ml_conical', 5)

    # Temperature modules
    temp_module_pcrtubes = protocol.load_module('temperature module gen2', 10)
    pcrtubes_cool = temp_module_pcrtubes.load_labware('opentrons_96_aluminumblock_generic_pcr_strip_200ul')
    temp_module_pcrtubes.set_temperature(4)

    temp_module_eppendorftubes = protocol.load_module('temperature module gen2', 7)
    eppendorftubes_cool = temp_module_eppendorftubes.load_labware('opentrons_24_aluminumblock_nest_2ml_snapcap')
    temp_module_eppendorftubes.set_temperature(4)  
    
    ## Define start reagents

    # In Eppendorf module
    DNA = eppendorftubes_cool.wells_by_name()["B6"]          # X mL 50 ng/muL DNA
    BufferW = eppendorftubes_cool.wells_by_name()["C5"]      # X mL
    Lysate = eppendorftubes_cool.wells_by_name()["C6"]       # X mL
    MQ = eppendorftubes_cool.wells_by_name()["D6"]           # 2 mL
    
    # In PCR module
    Mg_glut = pcrtubes_cool.wells_by_name()["A8"]      # X muL 1M
    K_glut = pcrtubes_cool.wells_by_name()["B8"]       # X muL 2M
    PEG8000 = pcrtubes_cool.wells_by_name()["C8"]      # muL 40%

    # In rack
    MM = rack.wells_by_name()["A1"]

    def find_index(lst,item):
        if item in lst:
            return lst.index(item)
        else:
            return None   

    def sorted_factors_to_lst(factors,reagent):
        """Extract factor values from pandas dataframe.
        Input is a string with name of reagent and outputs a sorted list"""
        lst = factors[reagent].tolist()
        lst.sort()
        return lst

    def pipette_vicious(mode):
        """Decrease flow rate for pipettes suitable for vicous liquids."""
        if mode == 0:     # OFF
                p20.flow_rate.aspirate = 7.56
                p20.flow_rate.dispense = 7.56
                p300.flow_rate.aspirate = 92.86
                p300.flow_rate.dispense = 92.86
        if mode == 1:     # ON
                p20.flow_rate.aspirate = 2
                p20.flow_rate.dispense = 2
                p300.flow_rate.aspirate = 10
                p300.flow_rate.dispense = 10

    def load_control(factors):
        """Load internal control wo/ DNA to well plate 
        w/ 3 mM Mg-glutamate, 60 mM K-glutamate, and 2% PEG-8000"""

        # No DNA
        p20.transfer(0.5, MQ, plate.wells()[nsamples], touch_tip=True)

        # Add reference factors
        Mg_final_conc = sorted_factors_to_lst(factors,'Mg-glutamate')
        K_final_conc = sorted_factors_to_lst(factors,'K-glutamate')
        P_final_conc = sorted_factors_to_lst(factors,'PEG-8000')
        if find_index(Mg_final_conc,3) is not None:
            dil_well = "A{}".format(find_index(Mg_final_conc,3)+1)
            p20.transfer(0.5, pcrtubes_cool.wells_by_name()[dil_well], plate.wells()[nsamples], touch_tip=True)
        else:
            # Make dilution and add to well
            p20.transfer(0.6, Mg_glut, pcrtubes_cool.wells_by_name()['C1'], touch_tip=True)
            p20.transfer(9.4, MQ, pcrtubes_cool.wells_by_name()['C1'], touch_tip=True)
            p20.transfer(0.5, pcrtubes_cool.wells_by_name()['C1'], plate.wells()[nsamples], touch_tip=True)
        if find_index(K_final_conc,60) is not None:
            dil_well = "B{}".format(find_index(K_final_conc,60)+1)
            p20.transfer(0.75, pcrtubes_cool.wells_by_name()[dil_well], plate.wells()[nsamples], touch_tip=True)
        else:
            # Make dilution and add to well
            p20.transfer(4, K_glut, pcrtubes_cool.wells_by_name()['C2'], touch_tip=True)
            p20.transfer(6, MQ, pcrtubes_cool.wells_by_name()['C2'], touch_tip=True)
            p20.transfer(0.75, pcrtubes_cool.wells_by_name()['C2'], plate.wells()[nsamples], touch_tip=True)
        if find_index(P_final_conc,2) is not None:
            dil_well = "C{}".format(find_index(P_final_conc,2)+1)
            p20.transfer(1.5, pcrtubes_cool.wells_by_name()[dil_well], plate.wells()[nsamples], touch_tip=True)
        else:
            # Make dilution and add to well
            p20.transfer(3.3, PEG8000, pcrtubes_cool.wells_by_name()['C3'], touch_tip=True)
            p20.transfer(6.6, MQ, pcrtubes_cool.wells_by_name()['C3'], touch_tip=True)
            p20.transfer(1.5, pcrtubes_cool.wells_by_name()['C3'], plate.wells()[nsamples], touch_tip=True)
        
    def calc_volume(reagent, final_conc_lst, stock_conc, final_vol, dil_factor):
        """
        Calculate volumes of reagent stock/dilution to get final concentration.
        Takes in a list of final concentrations, final volume, name of reagent and dilution factor
        and outputs two lists with volumes of stock/dilution and MilliQ.
        """
        reagent_vol_lst = []
        mq_vol_lst = []
        for i in range(len(final_conc_lst)):
            vol = ((final_conc_lst[i]*dil_factor)*final_vol)/stock_conc
            if vol > final_vol:
                raise ValueError("Error: Stock volume exceeds 30 muL. Change {} stock concentration.".format(reagent))
            elif vol < 0.5 and vol != 0:
                raise ValueError("Error: Stock volume lower than minimum required volume 0.5 muL. Change {} stock concentration.".format(reagent))
            else:
                if final_vol == ceil(vol):
                    reagent_vol_lst.append(ceil(vol))
                    mq_vol_lst.append(0)
                else:
                    if final_vol-vol < 0.5:
                        raise ValueError("Error: MilliQ volume lower than minimum required volume 0.5 muL. Change {} stock concentration.".format(reagent))
                    else:
                        reagent_vol_lst.append(vol)
                        mq_vol_lst.append(final_vol-vol)
        return reagent_vol_lst, mq_vol_lst

    def factors_dilution(factors):
        """Make dilutions of Mg-glut, Mg-glut and PEG-8000
        based on csv file input."""
        
        # Load reagent's stock solution concentration
        stock_conc = pd.read_csv('stocks.csv', header=0)
        Mg_stock_conc = stock_conc['Mg-glutamate (mM)'].iloc[0]
        K_stock_conc = stock_conc['K-glutamate (mM)'].iloc[0]
        P_stock_conc = stock_conc['PEG-8000 (%)'].iloc[0]

        # Load reagent's final concentration
        Mg_final_conc = sorted_factors_to_lst(factors,'Mg-glutamate')
        K_final_conc = sorted_factors_to_lst(factors,'K-glutamate')
        P_final_conc = sorted_factors_to_lst(factors,'PEG-8000')

        # Prepare dilutions
        for i in range(3):
            row = pcrtubes_cool.rows()[i]
            if i == 0:    # Mg-glutamate
                (Mg_glut_vol_lst, MQ_vol_lst) = calc_volume("Mg-glutamate",Mg_final_conc,Mg_stock_conc,30,20) # final volume is 30 muL
                if any(vol > 20 for vol in MQ_vol_lst):
                    for j in range(len(MQ_vol_lst)):
                        if MQ_vol_lst[j] <= 20:
                            p20.transfer(MQ_vol_lst[j], MQ, row[j])
                        else:
                            p300.transfer(MQ_vol_lst[j], MQ, row[j])
                else:
                    p20.transfer(MQ_vol_lst, MQ, row[:len(Mg_final_conc)])
                if any(vol > 20 for vol in Mg_glut_vol_lst):
                    for j in range(len(Mg_glut_vol_lst)):
                        if Mg_glut_vol_lst[j] <= 20:
                            p20.transfer(Mg_glut_vol_lst[j], Mg_glut, row[j], touch_tip=True, mix_after=(5,15))
                        else:
                            p300.transfer(Mg_glut_vol_lst[j], Mg_glut, row[j], touch_tip=True, mix_after=(5,15))
                else:
                    p20.transfer(Mg_glut_vol_lst, Mg_glut, row[:len(Mg_final_conc)], touch_tip=True, mix_after=(5,15))


            elif i == 1:  # K-glutamate
                (K_glut_vol_lst, MQ_vol_lst) = calc_volume("K-glutamate",K_final_conc,K_stock_conc,30,13) # final volume is 30 muL
                if any(vol > 20 for vol in MQ_vol_lst):
                    for j in range(len(MQ_vol_lst)):
                        if MQ_vol_lst[j] <= 20:
                            p20.transfer(MQ_vol_lst[j], MQ, row[j])
                        else:
                            p300.transfer(MQ_vol_lst[j], MQ, row[j])
                else:
                    p20.transfer(MQ_vol_lst, MQ, row[:len(K_final_conc)])
                if any(vol > 20 for vol in K_glut_vol_lst):
                    for j in range(len(K_glut_vol_lst)):
                        if K_glut_vol_lst[j] <= 20:
                            p20.transfer(K_glut_vol_lst[j], K_glut, row[j], touch_tip=True, mix_after=(5,15))
                        else:
                            p300.transfer(K_glut_vol_lst[j], K_glut, row[j], touch_tip=True, mix_after=(5,15))
                else:
                    p20.transfer(K_glut_vol_lst, K_glut, row[:len(K_final_conc)], touch_tip=True, mix_after=(5,15))
            
            elif i == 2:  # PEG-8000
                (PEG8000_vol_lst, MQ_vol_lst) = calc_volume("PEG-8000",P_final_conc,P_stock_conc,30,6.66) # final volume is 30 muL
                if any(vol > 20 for vol in MQ_vol_lst):
                    for j in range(len(MQ_vol_lst)):
                        if MQ_vol_lst[j] <= 20:
                            p20.transfer(MQ_vol_lst[j], MQ, row[j])
                        else:
                            p300.transfer(MQ_vol_lst[j], MQ, row[j])
                else:
                    p20.transfer(MQ_vol_lst, MQ, row[:len(P_final_conc)])
                pipette_vicious(1)  # ON
                if any(vol > 20 for vol in PEG8000_vol_lst):
                    for j in range(len(PEG8000_vol_lst)):
                        if PEG8000_vol_lst[j] <= 20:
                            p20.transfer(PEG8000_vol_lst[j], PEG8000, row[j], touch_tip=True, mix_after=(5,15))
                        else:
                            p300.transfer(PEG8000_vol_lst[j], PEG8000, row[j], touch_tip=True, mix_after=(5,15))
                else:
                    p20.transfer(PEG8000_vol_lst, PEG8000, row[:len(P_final_conc)], touch_tip=True, mix_after=(5,15))
                pipette_vicious(0)  # OFF

    def cfe_mastermix_prep():
        """Prepare mastermix excl. factors to optimize.
        Mastermix consits of BufferW and lysate"""
         
        lysate_vol = 4 * (nsamples * 1.1)  #uL
        bufferW_vol = 3 * (nsamples * 1.1)  #uL
        p300.transfer(lysate_vol, Lysate, MM,  touch_tip=True, blow_out=True, blowout_location='source well')
        p300.transfer(bufferW_vol, BufferW, MM, blow_out=True, blowout_location='source well')

    def mix_mastermix():
        p300.pick_up_tip()
        for i in range(5):
            p300.aspirate(300, MM.bottom(1))
            p300.dispense(300, MM.top(-2))
        p300.drop_tip()    
    
    def load_combinations(factors):
        """Makes a full factorial experimental design for the 3 factors
        Mg-glutamate, K-glutamate and PEG-8000. This is used to load
        all possible combinations. All wells with the same condition 
        for a factor is loaded simultanously to the 384-well plate. 
        Input is the csv file with factors."""

        # Build full factorial design and sort randomly
        ff = build.full_fact(factors).sample(frac=1)
        ff.to_csv('DOE.csv')
        ff_indexed = ff.index
        Mg_final_conc = sorted_factors_to_lst(factors,'Mg-glutamate')
        K_final_conc = sorted_factors_to_lst(factors,'K-glutamate')
        P_final_conc = sorted_factors_to_lst(factors,'PEG-8000')

        # Add Mg-glutamate combinations
        for i in range(len(Mg_final_conc)):
            conc_i = ff["Mg-glutamate"]==Mg_final_conc[i]
            conc_i_wells = []
            [conc_i_wells.append(j) for j in range(len(conc_i)) if conc_i.tolist()[j] ]
            dil_well = "A{}".format(i+1)
            p20.distribute(0.5, pcrtubes_cool.wells_by_name()[dil_well], [plate.wells()[well] for well in conc_i_wells], touch_tip=True, blow_out=True, blowout_location='source well')

        # Add K-glutamate combinations
        for i in range(len(K_final_conc)):
            conc_i = ff["K-glutamate"]==K_final_conc[i]
            conc_i_wells = []
            [conc_i_wells.append(j) for j in range(len(conc_i)) if conc_i.tolist()[j] ]
            dil_well = "B{}".format(i+1)
            p20.distribute(0.75, pcrtubes_cool.wells_by_name()[dil_well], [plate.wells()[well] for well in conc_i_wells], touch_tip=True, blow_out=True, blowout_location='source well')
        
        # Add PEG-8000 combinations
        for i in range(len(P_final_conc)):
            conc_i = ff["PEG-8000"]==P_final_conc[i]
            conc_i_wells = []
            [conc_i_wells.append(j) for j in range(len(conc_i)) if conc_i.tolist()[j] ]
            dil_well = "C{}".format(i+1)
            p20.distribute(1.5, pcrtubes_cool.wells_by_name()[dil_well], [plate.wells()[well] for well in conc_i_wells], touch_tip=True, blow_out=True, blowout_location='source well')            


    ## Protocol workflow
    
    # Dilute factors based on data from csv file
    factors = pd.read_csv('factors.csv', header=0)
    factors_dilution(factors)
    
    
    # Prepare and load buffer mix excl. factors for optimization
    cfe_mastermix_prep()
    mix_mastermix()
    p300.distribute(7, MM, plate.wells()[:nsamples+1], touch_tip=True, blow_out=True, blowout_location='source well')

    # Load control wo/ DNA and reference concentration for all factors
    load_control(factors)

    # Load combinations of factors
    load_combinations(factors)

    # Add DNA to initiate cell-free expression
    # change pipette tip?
    p20.distribute(0.5, DNA, plate.wells()[:nsamples], touch_tip=True, mix_before=(3,15), blow_out=True, blowout_location='source well')

    temp_module_pcrtubes.deactivate()
    temp_module_eppendorftubes.deactivate()
    
