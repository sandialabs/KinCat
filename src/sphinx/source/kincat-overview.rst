KinCat Overview
===============

A basic KMC algorithm is below. This is often called the resident time algorithm (RTA) and is inherently sequential. In KinCat, we consider three solver implementations: 1) serial reference implementation, 2) sublattice parallel implementation, and 3) batch parallel implementation.

.. code-block:: c++

   SolverSerialRTA::advance(in: t, in: t_step, in: n_max_iterations,
                            in: lattice, in: dictionary, in: rates, in: rates_scan) { 		
     /// update event rates on lattice 
     updateEventRatesLattice(in: lattice, in: dictionary, out: rates);
     
     for (iter=0; iter<n_max_iteration && t<t_end; ++iter) {
       /// perform prefix sum
       scanRatesLattice(in: rates, out: rates_scan);

       getRandomNumber(out: pi, out: zeta);
       sum_rates = rates_scan(in: rates.extent(0));
       dt = log(zeta)/sum_rates;
       rate_to_search = sum_rates * pi;
  
       /// cid is the cell index corresponding to a specific rate (sum_rates*pi) in the scanned array 
       searchCellIndex(in: rate_to_search, in: rates_scan, out: cid);

       /// select an event in the selected cell 
       findEvent(in: cid, in: dictionary, out: event);

       /// update the cell configuration
       updateLattice(in: cid, in: event, out: lattice);

       /// local updates of rates on the neighborhood of the selected cell
       updateEventRates(in: lattice, in: dictionary, out: rates);
  
       /// advance time
       t += dt;
     }
   }

The sublattice algorithm decomposes the lattice into multiple subdomains. In each subdomain, a sequential RTA algorithm runs and time is updated asynchronously. To avoid the potential conflicts on the boundary regions, the algorithm updates even and odd numbered domains separately.

.. code-block:: c++

   SolverSublattice::advance(in: t, in: t_step, in: n_max_iterations,
                             in: lattice, in: dictionary, in: rates, in: rates_scan) { 		
		
   /// update event rates on lattice 
   updateEventRatesLattice(in: lattice, in: dictionary, out: rates);

   /// even odd selector
   quad = { {0,0}, {1,0), {0,1}, {1,1} };
   
   for (iter=0; iter<n_max_iteration && t<t_end; ++iter) {
     for (q=0; q<4; ++q) {
       parallel_for(in: team_policy(n_domains), in: LAMBDA(member) {
         /// compute domain indices corresponding to the member's league rank
         getDomainIndex(in: member, in: lattice, out: d0, out: d1);
	 
	 /// select even or odd domain in 2D
	 if (d0%2 == quad[q][0] && d1%2 == quad[q][1] && t(d0,d1) < t_step) {
	   /// get domain specific containers
	   auto rates = rantes_lattice(d0, d1);	   
	   auto rates_scan = rantes_scan_lattice(d0, d1);

	   /// nested parallel scan over rates
	   parallelScanRates(in: member, in: rates, out: rates_scan);

	   /// update event 
	   single([=]() {
             getRandomNumber(out: pi, out: zeta);
             sum_rates = rates_scan(in: rates.extent(0));
             dt = log(zeta)/sum_rates;
             rate_to_search = sum_rates * pi;
            
             /// cid is the cell index corresponding to a specific rate (sum_rates*pi) in the scanned array 
             searchCellIndex(in: rate_to_search, in: rates_scan, out: cid);

             /// select an event in the selected cell 	     
             findEvent(in: cid, out: event);

             /// update the cell configuration	   
             updateLattice(in: cid, in: event, out: lattice);

	     /// update the domain clock
	     t(d0,d1) += dt;
	   });

           /// local updates of rates on the neighborhood of the selected cell	   
           parallelUpdateEventRates(in: member, in: lattice, in: dictionary, in: interaction_range, out: rates);
	 }
       });
     }
   }

The batch RTA solver runs the serial RTA, but in multiple independent simulations or samples. This may be desired for statistical or UQ analysis. The user interface and container structure needs to be changed to handle the use case correctly. For the batch mode, we implement a separate executable. 

.. autosummary::
   :toctree: generated
	     
