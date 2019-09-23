converter = AvgCnvFixedDC(690.0, 2750000.0) #S_rated goes in Watts
dc_source = FixedDCSource(600.0) #Not in the original data, guessed.
filter = LCLFilter(0.08, 0.003, 0.074, 0.2, 0.01)
pll = PLL(500.0, 0.084, 4.69)

virtual_H = VirtualInertia(2.0, 400.0, 20.0, 2*pi*50.0)
Q_control = ReactivePowerDroop(0.2, 1000.0)
outer_control = VirtualInertiaQdroop(virtual_H, Q_control)

vsc = CurrentControl(1.27, 14.3, 0.59, 736.0, 50.0, 0.5, 0.2,
0.00)

DynInverter(1,
            :DARCO,
            nodes_OMIB[2],
            1.0,
            1.0,
            0.4,
            0.0,
            2.75e3,
            converter,
            outer_control,
            vsc,
            dc_source,
            pll,
            filter)
