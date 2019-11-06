#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Top Block
# Generated: Mon Jan 28 22:15:23 2019
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from PyQt4 import Qt
from gnuradio import blocks
from gnuradio import channels
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import qtgui
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import pmt
import sip
import sys
from gnuradio import qtgui


class top_block(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Top Block")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Top Block")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "top_block")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())


        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 32000
        self.noise_amp = noise_amp = 0

        ##################################################
        # Blocks
        ##################################################
        self.qtgui_sink_x_0_0_0 = qtgui.sink_c(
        	1024, #fftsize
        	firdes.WIN_BLACKMAN_hARRIS, #wintype
        	0, #fc
        	samp_rate, #bw
        	"", #name
        	True, #plotfreq
        	True, #plotwaterfall
        	True, #plottime
        	True, #plotconst
        )
        self.qtgui_sink_x_0_0_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_0_0_win = sip.wrapinstance(self.qtgui_sink_x_0_0_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_sink_x_0_0_0_win)

        self.qtgui_sink_x_0_0_0.enable_rf_freq(False)



        self.qtgui_sink_x_0_0 = qtgui.sink_c(
        	1024, #fftsize
        	firdes.WIN_BLACKMAN_hARRIS, #wintype
        	0, #fc
        	samp_rate, #bw
        	"", #name
        	True, #plotfreq
        	True, #plotwaterfall
        	True, #plottime
        	True, #plotconst
        )
        self.qtgui_sink_x_0_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_0_win = sip.wrapinstance(self.qtgui_sink_x_0_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_sink_x_0_0_win)

        self.qtgui_sink_x_0_0.enable_rf_freq(False)



        self.qtgui_sink_x_0 = qtgui.sink_c(
        	1024, #fftsize
        	firdes.WIN_BLACKMAN_hARRIS, #wintype
        	0, #fc
        	samp_rate, #bw
        	"", #name
        	True, #plotfreq
        	True, #plotwaterfall
        	True, #plottime
        	True, #plotconst
        )
        self.qtgui_sink_x_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_win = sip.wrapinstance(self.qtgui_sink_x_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_sink_x_0_win)

        self.qtgui_sink_x_0.enable_rf_freq(False)



        self.channels_channel_model_0_1_0 = channels.channel_model(
        	noise_voltage=noise_amp,
        	frequency_offset=0.0,
        	epsilon=1.0,
        	taps=(0.1- 0.1j, ),
        	noise_seed=0,
        	block_tags=False
        )
        self.channels_channel_model_0_1 = channels.channel_model(
        	noise_voltage=noise_amp,
        	frequency_offset=0.0,
        	epsilon=1.0,
        	taps=(1- 1j, ),
        	noise_seed=0,
        	block_tags=False
        )
        self.channels_channel_model_0_0_0_0 = channels.channel_model(
        	noise_voltage=noise_amp,
        	frequency_offset=0.0,
        	epsilon=1.0,
        	taps=(0.2-0.2j, ),
        	noise_seed=0,
        	block_tags=False
        )
        self.channels_channel_model_0_0_0 = channels.channel_model(
        	noise_voltage=noise_amp,
        	frequency_offset=0.0,
        	epsilon=1.0,
        	taps=(1.1-1.1j, ),
        	noise_seed=0,
        	block_tags=False
        )
        self.channels_channel_model_0_0 = channels.channel_model(
        	noise_voltage=noise_amp,
        	frequency_offset=0.0,
        	epsilon=1.0,
        	taps=(1, ),
        	noise_seed=0,
        	block_tags=False
        )
        self.channels_channel_model_0 = channels.channel_model(
        	noise_voltage=noise_amp,
        	frequency_offset=0.0,
        	epsilon=1.0,
        	taps=(5, ),
        	noise_seed=0,
        	block_tags=False
        )
        self.blocks_file_source_1 = blocks.file_source(gr.sizeof_gr_complex*1, 'C:\\Users\\pedra\\Desktop\\NOMA_2X3\\signals\\super_tx_2.dat', True)
        self.blocks_file_source_1.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_gr_complex*1, 'C:\\Users\\pedra\\Desktop\\NOMA_2X3\\signals\\super_tx_1.dat', True)
        self.blocks_file_source_0.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_sink_0_0_0 = blocks.file_sink(gr.sizeof_gr_complex*1, 'C:\\Users\\pedra\\Desktop\\NOMA_2X3\\signals\\noma_rx_u1.dat', False)
        self.blocks_file_sink_0_0_0.set_unbuffered(False)
        self.blocks_file_sink_0_0 = blocks.file_sink(gr.sizeof_gr_complex*1, 'C:\\Users\\pedra\\Desktop\\NOMA_2X3\\signals\\noma_rx_u2.dat', False)
        self.blocks_file_sink_0_0.set_unbuffered(False)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_gr_complex*1, 'C:\\Users\\pedra\\Desktop\\NOMA_2X3\\signals\\noma_rx_u3.dat', False)
        self.blocks_file_sink_0.set_unbuffered(False)
        self.blocks_add_xx_0_0_0 = blocks.add_vcc(1)
        self.blocks_add_xx_0_0 = blocks.add_vcc(1)
        self.blocks_add_xx_0 = blocks.add_vcc(1)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_add_xx_0, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.blocks_add_xx_0, 0), (self.qtgui_sink_x_0, 0))
        self.connect((self.blocks_add_xx_0_0, 0), (self.blocks_file_sink_0_0, 0))
        self.connect((self.blocks_add_xx_0_0, 0), (self.qtgui_sink_x_0_0, 0))
        self.connect((self.blocks_add_xx_0_0_0, 0), (self.blocks_file_sink_0_0_0, 0))
        self.connect((self.blocks_add_xx_0_0_0, 0), (self.qtgui_sink_x_0_0_0, 0))
        self.connect((self.blocks_file_source_0, 0), (self.channels_channel_model_0, 0))
        self.connect((self.blocks_file_source_0, 0), (self.channels_channel_model_0_1, 0))
        self.connect((self.blocks_file_source_0, 0), (self.channels_channel_model_0_1_0, 0))
        self.connect((self.blocks_file_source_1, 0), (self.channels_channel_model_0_0, 0))
        self.connect((self.blocks_file_source_1, 0), (self.channels_channel_model_0_0_0, 0))
        self.connect((self.blocks_file_source_1, 0), (self.channels_channel_model_0_0_0_0, 0))
        self.connect((self.channels_channel_model_0, 0), (self.blocks_add_xx_0, 0))
        self.connect((self.channels_channel_model_0_0, 0), (self.blocks_add_xx_0, 1))
        self.connect((self.channels_channel_model_0_0_0, 0), (self.blocks_add_xx_0_0, 1))
        self.connect((self.channels_channel_model_0_0_0_0, 0), (self.blocks_add_xx_0_0_0, 1))
        self.connect((self.channels_channel_model_0_1, 0), (self.blocks_add_xx_0_0, 0))
        self.connect((self.channels_channel_model_0_1_0, 0), (self.blocks_add_xx_0_0_0, 0))

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "top_block")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.qtgui_sink_x_0_0_0.set_frequency_range(0, self.samp_rate)
        self.qtgui_sink_x_0_0.set_frequency_range(0, self.samp_rate)
        self.qtgui_sink_x_0.set_frequency_range(0, self.samp_rate)

    def get_noise_amp(self):
        return self.noise_amp

    def set_noise_amp(self, noise_amp):
        self.noise_amp = noise_amp
        self.channels_channel_model_0_1_0.set_noise_voltage(self.noise_amp)
        self.channels_channel_model_0_1.set_noise_voltage(self.noise_amp)
        self.channels_channel_model_0_0_0_0.set_noise_voltage(self.noise_amp)
        self.channels_channel_model_0_0_0.set_noise_voltage(self.noise_amp)
        self.channels_channel_model_0_0.set_noise_voltage(self.noise_amp)
        self.channels_channel_model_0.set_noise_voltage(self.noise_amp)


def main(top_block_cls=top_block, options=None):

    from distutils.version import StrictVersion
    if StrictVersion(Qt.qVersion()) >= StrictVersion("4.5.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()
    tb.start()
    tb.show()

    def quitting():
        tb.stop()
        tb.wait()
    qapp.connect(qapp, Qt.SIGNAL("aboutToQuit()"), quitting)
    qapp.exec_()


if __name__ == '__main__':
    main()
