# Model simulating the CO2, NO3-, and pH kinetics in a pH-stat based on empirical values, physical parameters in the reactor, 
# and cell-specific parameters including maximum growth rate and rate of electron transfer.

# version 2
# Marte Mølsæter Maråk and Lars Bakken, 23.07.24

# Import
import matplotlib.pyplot as plt
import pandas as pd
import math


def run_ph_stat_model(temp, od, co2_initial, ph_initial, no3_initial, n2_sparging, v_initial, k_feed, ph_setpoint,
                      f_fraction, time_step, total_duration):
    # Define experimental settings
    hno3_reservoir = 5  # M concentration in the reservoir
    hno3_injection_value = 0.23  # mL injected from the reservoir

    # Define organism specific parameters
    my_max = 0.198  # h-1
    y_n = 0.75  # mol cell-C mol-1 NO3 respired, Yield coefficient for cell carbon per mole of nitrate respired.
    cell_dw = 300  # fg DW cell-1, Cell dry weight in femtograms per cell.
    n_od = 1250000000  # cells mL-1 OD-1, Number of cells per mL per unit of optical density (OD).

    # Define dependent organism specific parameters
    y_e = y_n / 5  # mol cell-C mol-1 e-, Yield of cell carbon per mole of electrons.
    v_max = my_max / y_e  # mol e- mol-1 cell-C h-1, Maximum rate of electron transfer per mole of cell carbon per hour
    cell_c = cell_dw * 0.000000000000001 * 0.48 / 12  # mol C cell-1, 48% of DW is C, atomic mass of C = 12
    v_e_max = v_max * cell_c  # mol e- cell-1 h-1, Maximum rate of electron transfer per cell per hour.
    v_no3_max = v_e_max / 5  # mol NO3 cell-1 h-1,  Maximum rate of nitrate reduction per cell per hour.
    v_co2_no3 = 5 / 4  # mol CO2 mol-1 NO3 respired, Moles of CO2 produced per mole of NO3 respired
    km_no3 = 0.000002  # M NO2, Michaelis constant for nitrate in molarity (M).

    # Physical parameters
    kh_co2 = 0.03  # M atm-1, solubility of uncharged CO2 (Henrys constant)
    v_mol = 24.8  # L mol-1 at 30dC
    temp_k = temp + 273.15  # temperature in Kelvin
    cells_initial = od * n_od * 1000 * v_initial  # initial number of cells
    v_no3_max_reactor = cells_initial * v_no3_max / 3600  # M s-1 Vmax NO3 for the initial number of cells

    # Dissociation model - Appollo C.A, from KinCalc. mol L-1 at 1 atm partial pressure of CO2
    co2_model = 10 ** (108.3865 + 0.01985076 * temp_k - 6919.53 / temp_k - 40.45154 * math.log10(
        temp_k) + 669365 / (temp_k ** 2))

    # Empirical values
    d_co2 = 43  # deltapH M-1 CO2 total, pH depression by CO2, 43 at pH 7 and 70 at pH 7.4
    d_hno3 = 41  # deltapH M-1 NO3, pH depression by HNO3

    # Time axis
    time_axis = list(range(0, total_duration + time_step, time_step))

    # Initialize lists
    hno3_injection = [0]
    gluc_injection = [0]
    volume = [v_initial]
    no3_mol = [no3_initial * v_initial / 1000]  # Initial moles of NO3 based on initial concentration and volume
    no3_conc = [no3_mol[0] / volume[0]]  # Initial NO3 concentration
    v_no3 = [(v_no3_max_reactor * no3_mol[0]) / (no3_mol[0] + km_no3)]
    v_co2 = [v_no3[0] * v_co2_no3]
    co2_mol = [co2_initial]
    co2_total = [co2_mol[0] / volume[0]]
    p_h_change_due_to_co2 = [- d_co2 * co2_total[0]]
    p_h_change_due_to_hno3 = [- d_hno3 * no3_mol[0]]
    ph = [ph_initial + p_h_change_due_to_co2[0] + p_h_change_due_to_hno3[0]]
    hco3 = [co2_model * 10 ** (-356.3094 - 0.061 * temp_k + 21834.37 / temp_k + 126.8339 * math.log10(temp_k) - 1670000
                               / (temp_k ** 2)) / (10 ** -ph[0])]
    co32 = [hco3[0] * 10 ** (-107.8871 - 0.03252849 * temp_k + 5151.79 / temp_k + 38.92561 * math.log10(temp_k)
                             - 563713.9 / (temp_k ** 2)) / (10 ** -ph[0])]
    co2_total_model = [co2_model + hco3[0] + co32[0]]
    fraction_undissociated = [co2_model / co2_total_model[0]]
    fugacity = [co2_total[0] * fraction_undissociated[0] / kh_co2]  # Initialize Fugacity list
    t_co2 = [fugacity[0] * f_fraction * (n2_sparging / 60) / v_mol]  # Initialize T_CO2 list

    # initialize lists - cell calculations
    cells = [cells_initial]  # initialise list of cell number, cells
    cells_growth = [v_no3[0] * y_n / cell_c]  # initialise list of cell growth per time, cells -1
    growth_rate = [cells_growth[0] * 3600 / cells[0]]  # growth rate, h-1

    # Initialize lists - CO2 concentration in bursting bubbles
    co2_bubbles = [fugacity[0] * f_fraction * 100]  # Vol% CO2 in bubbles
    co2_headspace = [0]  # L CO2 in headspace
    co2_outflow = [100 * co2_headspace[0] / (4.3 - volume[0])]  # Vol% CO2 in outflow (total reactor volume 4.3 - liq v)

    # Iterate over time points
    for i in range(1, len(time_axis)):
        # Calculate pH-based HNO3 injection
        hno3_injection.append(hno3_injection_value / 1000 if ph[-1] > ph_setpoint else 0)

        # Total volume calculations
        gluc_injection.append(hno3_injection[-1] * k_feed)
        volume.append(volume[-1] + hno3_injection[-1] + gluc_injection[-1])

        # Forward Euler
        no3_mol_temp = no3_mol[-1] - time_step * v_no3[-1] + hno3_reservoir * hno3_injection[-1]
        no3_mol.append(max(no3_mol_temp, 0))  # Ensure NO3_mol does not go below 0
        co2_mol.append(co2_mol[-1] + time_step * (v_co2[-1] - t_co2[-1]))

        # NO3 and CO2 mol calculations
        no3_conc.append(no3_mol[-1] / volume[-1])
        v_no3.append((v_no3_max_reactor * no3_mol[-1]) / (no3_mol[-1] + km_no3))
        v_co2.append(v_no3[-1] * v_co2_no3)
        co2_total.append(co2_mol[-1] / volume[-1])

        # Calculate current pH
        p_h_change_due_to_co2.append(-d_co2 * co2_total[-1])
        p_h_change_due_to_hno3.append(-d_hno3 * no3_mol[-1])
        ph.append(ph_initial + p_h_change_due_to_co2[-1] + p_h_change_due_to_hno3[-1])

        # CO2 modeling
        hco3.append(co2_model * 10 ** (-356.3094 - 0.061 * temp_k + 21834.37 /
                                       temp_k + 126.8339 * math.log10(temp_k) - 1670000 / (temp_k ** 2)) / (
                            10 ** -ph[-1]))
        co32.append(hco3[-1] * 10 ** (-107.8871 - 0.03252849 * temp_k + 5151.79 / temp_k + 38.92561 *
                                      math.log10(temp_k) - 563713.9 / (temp_k ** 2)) / (10 ** -ph[-1]))
        co2_total_model.append(co2_model + hco3[-1] + co32[-1])
        fraction_undissociated.append(co2_model / co2_total_model[-1])

        # CO2 transport
        fugacity.append(co2_total[-1] * fraction_undissociated[-1] / kh_co2)
        t_co2.append(fugacity[-1] * f_fraction * (n2_sparging / 60) / v_mol)

        # Biomass calculations
        cells.append(cells[-1] + cells_growth[-1] * time_step)
        cells_growth.append(v_no3[-1] * y_n / cell_c)
        growth_rate.append(cells_growth[-1] * 3600 / cells[-1])

        # CO2 concentration in bursting bubbles and in the outflow
        co2_bubbles.append(fugacity[-1] * f_fraction * 100)
        co2_headspace.append(co2_headspace[-1] + (co2_bubbles[-1] / 100) * (n2_sparging / 60) * time_step -
                             co2_headspace[-1] * n2_sparging / (4.3 - volume[-2]) * time_step / 60)
        co2_outflow.append(100 * co2_headspace[-1] / (4.3 - volume[-1]))

    # Plotting the results
    time_axis_min = [t / 60 for t in time_axis]
    no3_conc_mM = [n * 1000 for n in no3_conc]
    co2_total_mM = [c * 1000 for c in co2_total]

    # Create a figure with subplots
    fig, axs = plt.subplots(2, 2, figsize=(18, 10))

    # Plot NO3_mol and CO2_mol on the left y-axis, and pH on the right y-axis
    ax1 = axs[0, 0]
    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('NO3_mol and CO2_mol (mol)')
    ax1.plot(time_axis_min, no3_conc_mM, label='NO3 [mM]', color='tab:red')
    ax1.plot(time_axis_min, co2_total_mM, label='CO2 [mM]', color='tab:blue')
    ax1.tick_params(axis='y')

    # Create a secondary y-axis for pH
    ax2 = ax1.twinx()
    ax2.set_ylabel('pH', color='orange')
    ax2.plot(time_axis_min, ph, label='pH', color='orange')
    ax2.tick_params(axis='y', labelcolor='orange')

    # Adding legends
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')

    ax1.set_title('Simulation')
    ax1.grid(True)

    # Plot the biomass calculations
    ax3 = axs[0, 1]
    ax4 = ax3.twinx()
    ax3.plot(time_axis_min, cells, 'g-', label='Cell Number')
    ax4.plot(time_axis_min, growth_rate, 'b-', label='Growth Rate')

    ax3.set_xlabel('Time (min)')
    ax3.set_ylabel('Cell Number', color='g')
    ax4.set_ylabel('Growth Rate (h-1)', color='b')

    ax3.legend(loc='upper left')
    ax4.legend(loc='upper right')
    ax3.set_title('Biomass Calculations')
    ax3.grid(True)

    # Plot the NO3 concentration over time
    axs[1, 0].plot(time_axis_min, no3_conc_mM, label='NO3 Concentration (mM)', color='tab:red')
    axs[1, 0].set_xlabel('Time (min)')
    axs[1, 0].set_ylabel('NO3 Concentration (mM)')
    axs[1, 0].set_title('NO3 Concentration')
    axs[1, 0].legend()
    axs[1, 0].grid(True)

    # Plot the CO2 calculations
    axs[1, 1].plot(time_axis_min, co2_bubbles, label='CO2 in Bubbles (Vol%)')
    axs[1, 1].plot(time_axis_min, co2_outflow, label='CO2 in Outflow (Vol%)')
    axs[1, 1].set_xlabel('Time (min)')
    axs[1, 1].set_ylabel('CO2 Concentration (Vol%)')
    axs[1, 1].set_title('CO2 Calculations')
    axs[1, 1].legend()
    axs[1, 1].grid(True)

    # Adjust layout to prevent overlap
    plt.tight_layout()
    plt.show()

    return time_axis_min, no3_conc_mM, co2_total_mM, ph, cells, growth_rate, co2_bubbles, co2_outflow


def export_to_excel(filename, time_axis_min, no3_conc_mM, co2_total_mM, ph, cells, growth_rate, co2_bubbles,
                    co2_outflow):
    # Create a DataFrame with the collected data
    data = {
        'Time (min)': time_axis_min,
        'NO3 Concentration (mM)': no3_conc_mM,
        'CO2 Concentration (mM)': co2_total_mM,
        'pH': ph,
        'Cell Number': cells,
        'Growth Rate (h-1)': growth_rate,
        'CO2 in Bubbles (Vol%)': co2_bubbles,
        'CO2 in Outflow (Vol%)': co2_outflow
    }

    df = pd.DataFrame(data)

    # Export DataFrame to Excel file
    df.to_excel(filename, index=False)
    print(f"Data exported to {filename}")


if __name__ == "__main__":
    # Run the model and get the results
    time_axis_min, no3_conc_mM, co2_total_mM, ph, cells, growth_rate, co2_bubbles, co2_outflow = run_ph_stat_model(
        temp=30, od=70, co2_initial=0.002, ph_initial=7.38, no3_initial=1, n2_sparging=0.5, v_initial=2,
        k_feed=0.6598, ph_setpoint=7, f_fraction=0.9, time_step=2, total_duration=60 * 10)
    # Export results to an Excel file
    export_to_excel('Test1_1.xlsx', time_axis_min, no3_conc_mM, co2_total_mM, ph, cells, growth_rate,
                    co2_bubbles, co2_outflow)
