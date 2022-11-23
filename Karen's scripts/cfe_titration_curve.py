from opentrons import protocol_api
import sys

metadata = {
    'apiLevel': '2.8',
    'protocolName': 'Cell Free expression titration curve',
    'description': '''This protocol is for setting up a 
     logaritmic serial dilution of a single reagent for cell free expression. 
     In total 7 concentrations of the reagnet are tested in technical triplicates
     including an internal control without DNA.''',
    'author': 'Karen Therkelsen (s173684@dtu.dk)',
}

## User inputs

# Define row in well-plate to load
row = "H"

# Define reagent stock concentration
conc = float(6) # mM

####################################################################################

# Display titration curve
conc_lst = []
for i in range(7):
    if i == 0:
        new_conc = float(conc)/2
        conc_lst.append(str(conc))
    new_conc = float(new_conc)/2
    conc_lst.append(str(new_conc))
print("The following titration row will be made:\n{}\n".format(",".join(conc_lst)))

nsamples = 8 * 3


def run(protocol):
    
    ## Load instrument, modules and labware ##

    #--------------#--------------#--------------#
    # PCR strips,  #   384-well   #    Trash     #
    # 4C   	   #   plate      #              #
    #--------------#--------------#--------------#
    # 2mLEppendorf #   P20 tips   #  P300 tips   #
    #  tubes, 4C   #              #              #
    #------------- #------------- #------------- #
    #              #              #              #
    #              #              #              #
    #------------- #------------- #------------- #
    #              #              #              #
    #              #              #              #
    #------------- #------------- #------------- #


    # Pipettes and tips
    tips20 = protocol.load_labware('opentrons_96_tiprack_20ul', 8)
    tips300 = protocol.load_labware('opentrons_96_tiprack_300ul', 9)
    p20 = protocol.load_instrument('p20_single_gen2', mount='right', tip_racks=[tips20])
    p300 = protocol.load_instrument('p300_single', mount='left', tip_racks=[tips300])

    # 384-well plate
    plate = protocol.load_labware('corning_384_wellplate_112ul_flat', 11)

    # Temperature modules
    temp_module_pcrtubes = protocol.load_module('temperature module gen2', 10)
    pcrtubes_cool = temp_module_pcrtubes.load_labware('opentrons_96_aluminumblock_generic_pcr_strip_200ul')

    temp_module_eppendorftubes = protocol.load_module('temperature module gen2', 4)
    eppendorftubes_cool = temp_module_eppendorftubes.load_labware('opentrons_24_aluminumblock_nest_2ml_snapcap')
    

    ## Define start reagents

    # In PCR module
    X = pcrtubes_cool.wells()[0]            # 15 uL
    DNA = pcrtubes_cool.wells()[8]          # 15 uL
    rNTP = pcrtubes_cool.wells()[9]         # 15 uL
    Lysate = pcrtubes_cool.wells()[10]      # 110 uL
    Buffer = pcrtubes_cool.wells()[11]      # 125 uL

    # In Eppendorf module
    MQ = eppendorftubes_cool.wells_by_name()["A1"]     # 1 mL
    MM = eppendorftubes_cool.wells_by_name()["A2"]


    ## Functions

    def serial_dilution():
        """Prepare logaritmic serial dilution with reagent"""
        p20.distribute(9, MQ, pcrtubes_cool.wells()[1:7])
        p20.pick_up_tip()
        p20.transfer(1, pcrtubes_cool.wells()[:6], pcrtubes_cool.wells()[1:7], mix_after=(3,5), new_tip="never")
        p20.drop_tip()

    def mix_mastermix():
        """Ensure homogenous mastermix before loding"""
        p300.pick_up_tip()
        for i in range(5):
            p300.aspirate(100, MM.bottom(1))
            p300.dispense(100, MM.top(-20))
        p300.drop_tip()
     
    def cfe_mastermix_prep():
        """Prepare CFE master mix excl. reagent"""
        lysate_vol = 4 * (nsamples * 1.2) #uL
        buffer_vol = 4.5 * (nsamples * 1.2)  #uL
        rNTP_vol = 0.5 * (nsamples * 1.2)  #uL

        p300.transfer(lysate_vol, Lysate, MM, touch_tip=True)
        p300.transfer(buffer_vol, Buffer, MM, touch_tip=True)
        p20.transfer(rNTP_vol, rNTP, MM, touch_tip=True)
    

    ## Protocol workflow
    
    #temp_module_pcrtubes.set_temperature(4)
    #temp_module_eppendorftubes.set_temperature(4)

    # Prepare master mix and serial dilution
    serial_dilution()
    cfe_mastermix_prep()

    # Add Add mastermix
    mix_mastermix()
    p20.distribute(9.0, MM, plate.rows_by_name()[row][:nsamples], touch_tip=True, blow_out=True, blowout_location='source well')
    
    # Add reagent
    for i in range(0,nsamples-3,3):
        dil_no = int((i+1)/3)
        if dil_no == 0:  # Load to internal control
            p20.pick_up_tip()
            p20.distribute(0.5, pcrtubes_cool.wells()[dil_no], plate.rows_by_name()[row][nsamples-3:nsamples], touch_tip=True, new_tip="never")
            p20.distribute(0.5, pcrtubes_cool.wells()[dil_no], plate.rows_by_name()[row][:i+3], touch_tip=True, new_tip="never")
            p20.drop_tip()
            p20.distribute(0.5, MQ, plate.rows_by_name()[row][nsamples-3:nsamples], touch_tip=True)
        else:
            p20.distribute(0.5, pcrtubes_cool.wells()[dil_no], plate.rows_by_name()[row][i:i+3], touch_tip=True)

    # Add DNA to initate CFE
    p20.distribute(0.5, DNA, plate.rows_by_name()[row][:nsamples-3], touch_tip=True)

    #temp_module_pcrtubes.deactivate()
    #temp_module_eppendorftubes.deactivate()
