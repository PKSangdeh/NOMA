# NOMA
This repo has our initial codes for offline implementation of NOMA with CVX

Purpose of the Test:

         o The test is for a MxN NOMA system where an M-antenna AP serves N single-antenna users.

         o The objective is to design the appropriate precoders at the AP to maximize the weighted
           sum-rate. This design should maintain fairness among users and also SIC requirements. 

         o This version of the project just is used as an initial test with offline processing and
           CVX solver for optimization. For real-time implementation, we used python and CVXOPT which 
           the codes will be uploaded after publishing the paper.

Procedure at high level:

         o Phase 1 (sounding): For the sounding, you can use both downlink or uplink channels. If 
           uplink is used, you should consider calibration as a must to compensate the DL/UL mismatch.
           The current version of the code uses DL channels directly to design the precoders.

         o Phase 2 (Precoder design): Upon obtaining the channel gains or vectors, the iterative 
           optimization algorithm calculates the precoders based on the channels and weights you
           define in the weighted sum-rate formula.

         o Phase 3 (Frame assembly): Based on the precoders, you generate a frame in which the payload
           part is the superimposition of all individual pre-coded payloads. The frame has an specific 
           format in the preamble part. Then the precoded frames are transmitted.

         o Phase 4 (SIC): The weakest user follows its regular manner, while others an iterative 
           interference cancellation method. There are two options for SIC, regular ZF-based SIC and
           our particular SIC which works slightly better.

Codes for each phase:

         o Phase 1: "tx_signal_gen" for generating sounding packets and "my_chan_estimate" for channel 
           estimation with user-defined over sampling.

         o Phase 2: CVX_PA/enum_t/find_tangs/double_check_sol/General_PA_main

         o Phase 3:  gen_super_frame or gen_super_frame_G_filter

         o Phase 4:  data_rx_decode or data_rx_decode_G_filter or ZF_vs_our_method based on the config.


Stored data: 

         o channel_ui.mat stores the downlink channels between the AP and user i at 52 subcarriers.

         o v_stack.mat has the stack of optimal precoders to be used for frame assembly.
         

Outputs:

         o The final results at users side are: EVM of decoded data, constellation of each SIC's stage,
           synchornization results, and freq. offsets.
