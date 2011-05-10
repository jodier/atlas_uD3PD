#############################################################################
# USER OPTIONS								    #
#############################################################################

AtlGeo = 'ATLAS-GEO-16-00-00'

CondDB = 'OFLCOND-SDR-BS7T-04-08'

#############################################################################

InputFormat = 'AOD'

isMC = True

#############################################################################

isEGamma = False

#############################################################################

InputFiles = [
	'../../AOD.280342._000152.pool.root'
]

OutputFile = 'output.root'

Stream = 'AANT'

#############################################################################
# JOB OPTIONS								    #
#############################################################################

from AthenaCommon.AppMgr import theApp
from AthenaCommon.AppMgr import ServiceMgr
from AthenaCommon.AppMgr import ToolSvc

#############################################################################

from AthenaCommon.GlobalFlags import globalflags

if isMC == False:
	globalflags.DataSource.set_Value_and_Lock('data')
else:
	globalflags.DataSource.set_Value_and_Lock('geant4')

globalflags.InputFormat.set_Value_and_Lock('pool')

globalflags.DetGeo.set_Value_and_Lock('atlas')

globalflags.DetDescrVersion.set_Value_and_Lock(AtlGeo)

globalflags.ConditionsTag.set_Value_and_Lock(CondDB)

#############################################################################

from AthenaCommon.AthenaCommonFlags import athenaCommonFlags

if InputFormat == 'ESD':
	athenaCommonFlags.PoolESDInput = InputFiles
if InputFormat == 'AOD':
	athenaCommonFlags.PoolAODInput = InputFiles

#############################################################################

import AthenaPoolCnvSvc.ReadAthenaPool

if InputFormat == 'ESD':
	ServiceMgr.EventSelector.InputCollections = athenaCommonFlags.PoolESDInput()
if InputFormat == 'AOD':
	ServiceMgr.EventSelector.InputCollections = athenaCommonFlags.PoolAODInput()

#############################################################################

from RecExConfig.RecFlags import rec

rec.doCBNT.set_Value_and_Lock(False)
rec.doWriteESD.set_Value_and_Lock(False)
rec.doWriteAOD.set_Value_and_Lock(False)
rec.doWriteTAG.set_Value_and_Lock(False)
rec.readESD.set_Value_and_Lock(InputFormat == 'ESD')
rec.readAOD.set_Value_and_Lock(InputFormat == 'AOD')
rec.doESD.set_Value_and_Lock(False)
rec.doAOD.set_Value_and_Lock(False)
rec.doDPD.set_Value_and_Lock(False)
rec.doHist.set_Value_and_Lock(False)

#############################################################################

from AthenaCommon.DetFlags import DetFlags

DetFlags.ID_setOn()
DetFlags.LAr_setOn()
DetFlags.Tile_setOn()
DetFlags.Calo_setOn()
DetFlags.Muon_setOn()

#############################################################################

from AthenaCommon.BFieldFlags import jobproperties

jobproperties.BField.solenoidOn.set_Value_and_Lock(True)
jobproperties.BField.barrelToroidOn.set_Value_and_Lock(True)
jobproperties.BField.endcapToroidOn.set_Value_and_Lock(True)

#############################################################################

from PoolSvc.PoolSvcConf import PoolSvc
ServiceMgr += PoolSvc(SortReplicas = True)

from DBReplicaSvc.DBReplicaSvcConf import DBReplicaSvc
ServiceMgr += DBReplicaSvc(UseCOOLSQLite = False)

#############################################################################

include('RecExCond/AllDet_detDescr.py')

#############################################################################

conddb.setGlobalTag(globalflags.ConditionsTag())

#############################################################################

from AtlasGeoModel import SetGeometryVersion

from AtlasGeoModel import GeoModelInit

#############################################################################

from AthenaCommon import CfgMgr

if not hasattr(ServiceMgr, 'THistSvc'):
	ServiceMgr += CfgMgr.THistSvc()

ServiceMgr.THistSvc.Output += ["%s DATAFILE='%s' TYP='ROOT' OPT='RECREATE'" % (Stream, OutputFile)]

#############################################################################
# TOOLS									    #
#############################################################################

from TriggerJobOpts.TriggerConfigGetter import TriggerConfigGetter

cfg = TriggerConfigGetter("ReadPool")

#################################
# TrigDecisionTool		#
#################################

from TrigDecisionTool.TrigDecisionToolConf import Trig__TrigDecisionTool

ToolSvc += Trig__TrigDecisionTool(name = 'TrigDecisionTool')

#################################
# AtlasExtrapolator		#
#################################

from TrkExTools.AtlasExtrapolator import AtlasExtrapolator

theAtlasExtrapolator = AtlasExtrapolator(name = 'AtlasExtrapolator')
ToolSvc += theAtlasExtrapolator

#################################
# ExtrapolateToCaloTool		#
#################################

from TrackToCalo.TrackToCaloConf import ExtrapolateToCaloTool

theExtrapolateToCaloTool = ExtrapolateToCaloTool(name = 'ExtrapolateToCaloTool', Extrapolator = theAtlasExtrapolator)
ToolSvc += theExtrapolateToCaloTool

#################################
# MCTruthClassifier		#
#################################

from MCTruthClassifier.MCTruthClassifierConf import MCTruthClassifier

ToolSvc += MCTruthClassifier(name = 'MCTruthClassifier', ExtrapolateToCaloTool = theExtrapolateToCaloTool)

#################################
# TrackToVertexIPEstimator	#
#################################

from TrkVertexFitterUtils.TrkVertexFitterUtilsConf import Trk__TrackToVertexIPEstimator

ToolSvc += Trk__TrackToVertexIPEstimator(name = 'TrackToVertexIPEstimator', Extrapolator = theAtlasExtrapolator)

#############################################################################
# ALGORITHM								    #
#############################################################################

import ROOT
import PyCintex

#############################################################################

import AthenaPython.PyAthena as PyAthena

#############################################################################

import array
import ctypes

#############################################################################

class uD3PD(PyAthena.Alg):

	#####################################################################

	def initialize(self):
		PyCintex.loadDict("libuD3PDDict")
		PyCintex.loadDict("libTrkTrackSummaryDict")

		self.StoreGateSvc = PyAthena.py_svc('StoreGateSvc')
		self.THistSvc = PyAthena.py_svc('THistSvc')


		if isEGamma:
			self.Tree1 = self.THistSvc['/%s/egamma'        % Stream] = ROOT.TTree('egamma'       , 'egamma'       )
			self.Tree2 = self.THistSvc['/%s/egammaTrigDec' % Stream] = ROOT.TTree('egammaTrigDec', 'egammaTrigDec')
		else:
			self.Tree1 = self.Tree2 = self.THistSvc['/%s/physics' % Stream] = ROOT.TTree('physics', 'physics')

		self.treeBuilder()
		self.treeCleaner()

		#########################
		# TOOLS			#
		#########################

		self.TrigDecisionTool = PyAthena.py_tool('Trig::TrigDecisionTool/TrigDecisionTool', iface = 'Trig::TrigDecisionTool')
		print('TrigDecisionTool:', self.TrigDecisionTool)

		if not self.TrigDecisionTool:
			return PyAthena.StatusCode.Failure

		##

		self.MCTruthClassifier = PyAthena.py_tool('MCTruthClassifier', iface = 'IMCTruthClassifier')
		print('MCTruthClassifier:', self.MCTruthClassifier)

		if not self.MCTruthClassifier:
			return PyAthena.StatusCode.Failure

		##

		self.TrackToVertexIPEstimator = PyAthena.py_tool('Trk::TrackToVertexIPEstimator/TrackToVertexIPEstimator', iface = 'Trk::ITrackToVertexIPEstimator')
		print('TrackToVertexIPEstimator:', self.TrackToVertexIPEstimator)

		if not self.TrackToVertexIPEstimator:
			return PyAthena.StatusCode.Failure

		##

		return PyAthena.StatusCode.Success

	#####################################################################

	def treeBuilder(self):
		#########################
		# EVENT			#
		#########################

		self.RunNumber = array.array('I', [0])
		self.EventNumber = array.array('I', [0])
		self.lbn = array.array('I', [0])

		##

		self.Tree1.Branch('RunNumber', self.RunNumber, 'RunNumber/i')
		self.Tree1.Branch('EventNumber', self.EventNumber, 'EventNumber/i')
		self.Tree1.Branch('lbn', self.lbn, 'lbn/i')

		if isEGamma:
			self.Tree2.Branch('RunNumber', self.RunNumber, 'RunNumber/i')
			self.Tree2.Branch('EventNumber', self.EventNumber, 'EventNumber/i')
			self.Tree2.Branch('lbn', self.lbn, 'lbn/i')

		#########################
		# PRIMARY VERTICES	#
		#########################

		self.vxp_n = array.array('i', [0])

		self.vxp_x = ROOT.std.vector(float)()
		self.vxp_y = ROOT.std.vector(float)()
		self.vxp_z = ROOT.std.vector(float)()

		self.vxp_nTracks = ROOT.std.vector(int)()

		##

		self.Tree1.Branch('vxp_n', self.vxp_n, 'vxp_n/I')

		self.Tree1.Branch('vxp_x', self.vxp_x)
		self.Tree1.Branch('vxp_y', self.vxp_y)
		self.Tree1.Branch('vxp_z', self.vxp_z)

		self.Tree1.Branch('vxp_nTracks', self.vxp_nTracks)

		#########################
		# ELECTRONS		#
		#########################

		self.el_n = array.array('i', [0])

		self.el_E = ROOT.std.vector(float)()
		self.el_Et = ROOT.std.vector(float)()
		self.el_pt = ROOT.std.vector(float)()
		self.el_eta = ROOT.std.vector(float)()
		self.el_phi = ROOT.std.vector(float)()
		self.el_charge = ROOT.std.vector(float)()
		self.el_author = ROOT.std.vector(int)()
		self.el_isEM = ROOT.std.vector(int)()
		self.el_OQ = ROOT.std.vector(int)()

		self.el_cl_E = ROOT.std.vector(float)()
		self.el_cl_pt = ROOT.std.vector(float)()
		self.el_cl_eta = ROOT.std.vector(float)()
		self.el_cl_phi = ROOT.std.vector(float)()

		self.el_Es2 = ROOT.std.vector(float)()
		self.el_etas2 = ROOT.std.vector(float)()
		self.el_phis2 = ROOT.std.vector(float)()

		self.el_loose = ROOT.std.vector(int)()
		self.el_medium = ROOT.std.vector(int)()
		self.el_tight = ROOT.std.vector(int)()

		self.el_reta = ROOT.std.vector(float)()
		self.el_weta2 = ROOT.std.vector(float)()

		self.el_ptcone20 = ROOT.std.vector(float)()
		self.el_ptcone30 = ROOT.std.vector(float)()
		self.el_ptcone40 = ROOT.std.vector(float)()

		self.el_Etcone20 = ROOT.std.vector(float)()
		self.el_Etcone30 = ROOT.std.vector(float)()
		self.el_Etcone40 = ROOT.std.vector(float)()

		self.el_Etcone20_pt_corrected = ROOT.std.vector(float)()
		self.el_Etcone30_pt_corrected = ROOT.std.vector(float)()
		self.el_Etcone40_pt_corrected = ROOT.std.vector(float)()

		self.el_nBLHits = ROOT.std.vector(int)()
		self.el_nPixHits = ROOT.std.vector(int)()
		self.el_nSCTHits = ROOT.std.vector(int)()
		self.el_nTRTHits = ROOT.std.vector(int)()

		self.el_trackd0 = ROOT.std.vector(float)()
		self.el_trackz0 = ROOT.std.vector(float)()
		self.el_trackphi = ROOT.std.vector(float)()
		self.el_tracktheta = ROOT.std.vector(float)()
		self.el_trackqoverp = ROOT.std.vector(float)()

		self.el_trackd0pvunbiased = ROOT.std.vector(float)()
		self.el_trackz0pvunbiased = ROOT.std.vector(float)()
		self.el_tracksigd0pvunbiased = ROOT.std.vector(float)()
		self.el_tracksigz0pvunbiased = ROOT.std.vector(float)()

		self.el_truth_type = ROOT.std.vector(int)()
		self.el_truth_mothertype = ROOT.std.vector(int)()
		self.el_truth_barbode = ROOT.std.vector(int)()
		self.el_truth_motherbarcode = ROOT.std.vector(int)()

		self.el_type = ROOT.std.vector(int)()
		self.el_origin = ROOT.std.vector(int)()
		self.el_typebkg = ROOT.std.vector(int)()
		self.el_originbkg = ROOT.std.vector(int)()

		##

		self.Tree1.Branch('el_n', self.el_n, 'el_n/I')

		self.Tree1.Branch('el_E', self.el_E)
		self.Tree1.Branch('el_Et', self.el_Et)
		self.Tree1.Branch('el_pt', self.el_pt)
		self.Tree1.Branch('el_eta', self.el_eta)
		self.Tree1.Branch('el_phi', self.el_phi)
		self.Tree1.Branch('el_charge', self.el_charge)
		self.Tree1.Branch('el_author', self.el_author)
		self.Tree1.Branch('el_isEM', self.el_isEM)
		self.Tree1.Branch('el_OQ', self.el_OQ)

		self.Tree1.Branch('el_cl_E', self.el_cl_E)
		self.Tree1.Branch('el_cl_pt', self.el_cl_pt)
		self.Tree1.Branch('el_cl_eta', self.el_cl_eta)
		self.Tree1.Branch('el_cl_phi', self.el_cl_phi)

		self.Tree1.Branch('el_Es2', self.el_Es2)
		self.Tree1.Branch('el_etas2', self.el_etas2)
		self.Tree1.Branch('el_phis2', self.el_phis2)

		self.Tree1.Branch('el_loose', self.el_loose)
		self.Tree1.Branch('el_medium', self.el_medium)
		self.Tree1.Branch('el_tight', self.el_tight)

		self.Tree1.Branch('el_reta', self.el_reta)
		self.Tree1.Branch('el_weta2', self.el_weta2)

		self.Tree1.Branch('el_ptcone20', self.el_ptcone20)
		self.Tree1.Branch('el_ptcone30', self.el_ptcone30)
		self.Tree1.Branch('el_ptcone40', self.el_ptcone40)

		self.Tree1.Branch('el_Etcone20', self.el_Etcone20)
		self.Tree1.Branch('el_Etcone30', self.el_Etcone30)
		self.Tree1.Branch('el_Etcone40', self.el_Etcone40)

		self.Tree1.Branch('el_Etcone20_pt_corrected', self.el_Etcone20_pt_corrected)
		self.Tree1.Branch('el_Etcone30_pt_corrected', self.el_Etcone30_pt_corrected)
		self.Tree1.Branch('el_Etcone40_pt_corrected', self.el_Etcone40_pt_corrected)

		self.Tree1.Branch('el_nBLHits', self.el_nBLHits)
		self.Tree1.Branch('el_nPixHits', self.el_nPixHits)
		self.Tree1.Branch('el_nSCTHits', self.el_nSCTHits)
		self.Tree1.Branch('el_nTRTHits', self.el_nTRTHits)

		self.Tree1.Branch('el_trackd0', self.el_trackd0)
		self.Tree1.Branch('el_trackz0', self.el_trackz0)
		self.Tree1.Branch('el_trackphi', self.el_trackphi)
		self.Tree1.Branch('el_tracktheta', self.el_tracktheta)
		self.Tree1.Branch('el_trackqoverp', self.el_trackqoverp)

		self.Tree1.Branch('el_trackd0pvunbiased', self.el_trackd0pvunbiased)
		self.Tree1.Branch('el_trackz0pvunbiased', self.el_trackz0pvunbiased)
		self.Tree1.Branch('el_tracksigd0pvunbiased', self.el_tracksigd0pvunbiased)
		self.Tree1.Branch('el_tracksigz0pvunbiased', self.el_tracksigz0pvunbiased)

		self.Tree1.Branch('el_truth_type', self.el_truth_type)
		self.Tree1.Branch('el_truth_mothertype', self.el_truth_mothertype)
		self.Tree1.Branch('el_truth_barbode', self.el_truth_barbode)
		self.Tree1.Branch('el_truth_motherbarcode', self.el_truth_motherbarcode)

		self.Tree1.Branch('el_type', self.el_type)
		self.Tree1.Branch('el_origin', self.el_origin)
		self.Tree1.Branch('el_typebkg', self.el_typebkg)
		self.Tree1.Branch('el_originbkg', self.el_originbkg)

		#########################
		# MUONS MUID		#
		#########################

		self.mu_muid_n = array.array('i', [0])

		self.mu_muid_E = ROOT.std.vector(float)()
		self.mu_muid_Et = ROOT.std.vector(float)()
		self.mu_muid_pt = ROOT.std.vector(float)()
		self.mu_muid_eta = ROOT.std.vector(float)()
		self.mu_muid_phi = ROOT.std.vector(float)()
		self.mu_muid_charge = ROOT.std.vector(float)()
		self.mu_muid_author = ROOT.std.vector(int)()

		self.mu_muid_loose = ROOT.std.vector(int)()
		self.mu_muid_medium = ROOT.std.vector(int)()
		self.mu_muid_tight = ROOT.std.vector(int)()

		self.mu_muid_ptcone20 = ROOT.std.vector(float)()
		self.mu_muid_ptcone30 = ROOT.std.vector(float)()
		self.mu_muid_ptcone40 = ROOT.std.vector(float)()

		self.mu_muid_etcone20 = ROOT.std.vector(float)()
		self.mu_muid_etcone30 = ROOT.std.vector(float)()
		self.mu_muid_etcone40 = ROOT.std.vector(float)()

		self.mu_muid_nBLHits = ROOT.std.vector(int)()
		self.mu_muid_nPixHits = ROOT.std.vector(int)()
		self.mu_muid_nSCTHits = ROOT.std.vector(int)()
		self.mu_muid_nTRTHits = ROOT.std.vector(int)()

		self.mu_muid_trackd0 = ROOT.std.vector(float)()
		self.mu_muid_trackz0 = ROOT.std.vector(float)()
		self.mu_muid_trackphi = ROOT.std.vector(float)()
		self.mu_muid_tracktheta = ROOT.std.vector(float)()
		self.mu_muid_trackqoverp = ROOT.std.vector(float)()

		self.mu_muid_trackd0pvunbiased = ROOT.std.vector(float)()
		self.mu_muid_trackz0pvunbiased = ROOT.std.vector(float)()
		self.mu_muid_tracksigd0pvunbiased = ROOT.std.vector(float)()
		self.mu_muid_tracksigz0pvunbiased = ROOT.std.vector(float)()

		self.mu_muid_truth_type = ROOT.std.vector(int)()
		self.mu_muid_truth_mothertype = ROOT.std.vector(int)()
		self.mu_muid_truth_barbode = ROOT.std.vector(int)()
		self.mu_muid_truth_motherbarcode = ROOT.std.vector(int)()

		##

		self.Tree1.Branch('mu_muid_n', self.mu_muid_n, 'mu_muid_n/I')

		self.Tree1.Branch('mu_muid_E', self.mu_muid_E)
		self.Tree1.Branch('mu_muid_Et', self.mu_muid_Et)
		self.Tree1.Branch('mu_muid_pt', self.mu_muid_pt)
		self.Tree1.Branch('mu_muid_eta', self.mu_muid_eta)
		self.Tree1.Branch('mu_muid_phi', self.mu_muid_phi)
		self.Tree1.Branch('mu_muid_charge', self.mu_muid_charge)
		self.Tree1.Branch('mu_muid_author', self.mu_muid_author)

		self.Tree1.Branch('mu_muid_loose', self.mu_muid_loose)
		self.Tree1.Branch('mu_muid_medium', self.mu_muid_medium)
		self.Tree1.Branch('mu_muid_tight', self.mu_muid_tight)

		self.Tree1.Branch('mu_muid_ptcone20', self.mu_muid_ptcone20)
		self.Tree1.Branch('mu_muid_ptcone30', self.mu_muid_ptcone30)
		self.Tree1.Branch('mu_muid_ptcone40', self.mu_muid_ptcone40)

		self.Tree1.Branch('mu_muid_etcone20', self.mu_muid_etcone20)
		self.Tree1.Branch('mu_muid_etcone30', self.mu_muid_etcone30)
		self.Tree1.Branch('mu_muid_etcone40', self.mu_muid_etcone40)

		self.Tree1.Branch('mu_muid_nBLHits', self.mu_muid_nBLHits)
		self.Tree1.Branch('mu_muid_nPixHits', self.mu_muid_nPixHits)
		self.Tree1.Branch('mu_muid_nSCTHits', self.mu_muid_nSCTHits)
		self.Tree1.Branch('mu_muid_nTRTHits', self.mu_muid_nTRTHits)

		self.Tree1.Branch('mu_muid_trackd0', self.mu_muid_trackd0)
		self.Tree1.Branch('mu_muid_trackz0', self.mu_muid_trackz0)
		self.Tree1.Branch('mu_muid_trackphi', self.mu_muid_trackphi)
		self.Tree1.Branch('mu_muid_tracktheta', self.mu_muid_tracktheta)
		self.Tree1.Branch('mu_muid_trackqoverp', self.mu_muid_trackqoverp)

		self.Tree1.Branch('mu_muid_trackd0pvunbiased', self.mu_muid_trackd0pvunbiased)
		self.Tree1.Branch('mu_muid_trackz0pvunbiased', self.mu_muid_trackz0pvunbiased)
		self.Tree1.Branch('mu_muid_tracksigd0pvunbiased', self.mu_muid_tracksigd0pvunbiased)
		self.Tree1.Branch('mu_muid_tracksigz0pvunbiased', self.mu_muid_tracksigz0pvunbiased)

		self.Tree1.Branch('mu_muid_truth_type', self.mu_muid_truth_type)
		self.Tree1.Branch('mu_muid_truth_mothertype', self.mu_muid_truth_mothertype)
		self.Tree1.Branch('mu_muid_truth_barbode', self.mu_muid_truth_barbode)
		self.Tree1.Branch('mu_muid_truth_motherbarcode', self.mu_muid_truth_motherbarcode)

		#########################
		# MUONS STACO		#
		#########################

		self.mu_staco_n = array.array('i', [0])

		self.mu_staco_E = ROOT.std.vector(float)()
		self.mu_staco_Et = ROOT.std.vector(float)()
		self.mu_staco_pt = ROOT.std.vector(float)()
		self.mu_staco_eta = ROOT.std.vector(float)()
		self.mu_staco_phi = ROOT.std.vector(float)()
		self.mu_staco_charge = ROOT.std.vector(float)()
		self.mu_staco_author = ROOT.std.vector(int)()

		self.mu_staco_loose = ROOT.std.vector(int)()
		self.mu_staco_medium = ROOT.std.vector(int)()
		self.mu_staco_tight = ROOT.std.vector(int)()

		self.mu_staco_ptcone20 = ROOT.std.vector(float)()
		self.mu_staco_ptcone30 = ROOT.std.vector(float)()
		self.mu_staco_ptcone40 = ROOT.std.vector(float)()

		self.mu_staco_etcone20 = ROOT.std.vector(float)()
		self.mu_staco_etcone30 = ROOT.std.vector(float)()
		self.mu_staco_etcone40 = ROOT.std.vector(float)()

		self.mu_staco_nBLHits = ROOT.std.vector(int)()
		self.mu_staco_nPixHits = ROOT.std.vector(int)()
		self.mu_staco_nSCTHits = ROOT.std.vector(int)()
		self.mu_staco_nTRTHits = ROOT.std.vector(int)()

		self.mu_staco_trackd0 = ROOT.std.vector(float)()
		self.mu_staco_trackz0 = ROOT.std.vector(float)()
		self.mu_staco_trackphi = ROOT.std.vector(float)()
		self.mu_staco_tracktheta = ROOT.std.vector(float)()
		self.mu_staco_trackqoverp = ROOT.std.vector(float)()

		self.mu_staco_trackd0pvunbiased = ROOT.std.vector(float)()
		self.mu_staco_trackz0pvunbiased = ROOT.std.vector(float)()
		self.mu_staco_tracksigd0pvunbiased = ROOT.std.vector(float)()
		self.mu_staco_tracksigz0pvunbiased = ROOT.std.vector(float)()

		self.mu_staco_truth_type = ROOT.std.vector(int)()
		self.mu_staco_truth_mothertype = ROOT.std.vector(int)()
		self.mu_staco_truth_barbode = ROOT.std.vector(int)()
		self.mu_staco_truth_motherbarcode = ROOT.std.vector(int)()

		##

		self.Tree1.Branch('mu_staco_n', self.mu_staco_n, 'mu_staco_n/I')

		self.Tree1.Branch('mu_staco_E', self.mu_staco_E)
		self.Tree1.Branch('mu_staco_Et', self.mu_staco_Et)
		self.Tree1.Branch('mu_staco_pt', self.mu_staco_pt)
		self.Tree1.Branch('mu_staco_eta', self.mu_staco_eta)
		self.Tree1.Branch('mu_staco_phi', self.mu_staco_phi)
		self.Tree1.Branch('mu_staco_charge', self.mu_staco_charge)
		self.Tree1.Branch('mu_staco_author', self.mu_staco_author)

		self.Tree1.Branch('mu_staco_loose', self.mu_staco_loose)
		self.Tree1.Branch('mu_staco_medium', self.mu_staco_medium)
		self.Tree1.Branch('mu_staco_tight', self.mu_staco_tight)

		self.Tree1.Branch('mu_staco_ptcone20', self.mu_staco_ptcone20)
		self.Tree1.Branch('mu_staco_ptcone30', self.mu_staco_ptcone30)
		self.Tree1.Branch('mu_staco_ptcone40', self.mu_staco_ptcone40)

		self.Tree1.Branch('mu_staco_etcone20', self.mu_staco_etcone20)
		self.Tree1.Branch('mu_staco_etcone30', self.mu_staco_etcone30)
		self.Tree1.Branch('mu_staco_etcone40', self.mu_staco_etcone40)

		self.Tree1.Branch('mu_staco_nBLHits', self.mu_staco_nBLHits)
		self.Tree1.Branch('mu_staco_nPixHits', self.mu_staco_nPixHits)
		self.Tree1.Branch('mu_staco_nSCTHits', self.mu_staco_nSCTHits)
		self.Tree1.Branch('mu_staco_nTRTHits', self.mu_staco_nTRTHits)

		self.Tree1.Branch('mu_staco_trackd0', self.mu_staco_trackd0)
		self.Tree1.Branch('mu_staco_trackz0', self.mu_staco_trackz0)
		self.Tree1.Branch('mu_staco_trackphi', self.mu_staco_trackphi)
		self.Tree1.Branch('mu_staco_tracktheta', self.mu_staco_tracktheta)
		self.Tree1.Branch('mu_staco_trackqoverp', self.mu_staco_trackqoverp)

		self.Tree1.Branch('mu_staco_trackd0pvunbiased', self.mu_staco_trackd0pvunbiased)
		self.Tree1.Branch('mu_staco_trackz0pvunbiased', self.mu_staco_trackz0pvunbiased)
		self.Tree1.Branch('mu_staco_tracksigd0pvunbiased', self.mu_staco_tracksigd0pvunbiased)
		self.Tree1.Branch('mu_staco_tracksigz0pvunbiased', self.mu_staco_tracksigz0pvunbiased)

		self.Tree1.Branch('mu_staco_truth_type', self.mu_staco_truth_type)
		self.Tree1.Branch('mu_staco_truth_mothertype', self.mu_staco_truth_mothertype)
		self.Tree1.Branch('mu_staco_truth_barbode', self.mu_staco_truth_barbode)
		self.Tree1.Branch('mu_staco_truth_motherbarcode', self.mu_staco_truth_motherbarcode)

		#########################
		# TRIGGERS		#
		#########################

		self.L1_EM14 = array.array('i', [0])
		self.L1_MU10 = array.array('i', [0])

		self.EF_e15_medium = array.array('i', [0])
		self.EF_e20_medium = array.array('i', [0])
		self.EF_e20_medium1 = array.array('i', [0])
		self.EF_2e12_medium = array.array('i', [0])

		self.EF_2g20_loose = array.array('i', [0])

		self.EF_mu10 = array.array('i', [0])
		self.EF_mu10_MG = array.array('i', [0])
		self.EF_mu13_MG = array.array('i', [0])
		self.EF_mu13_MG_tight = array.array('i', [0])
		self.EF_mu20_MG = array.array('i', [0])

		##

		self.Tree2.Branch('L1_EM14', self.L1_EM14, 'L1_EM14/I')
		self.Tree2.Branch('L1_MU10', self.L1_MU10, 'L1_MU10/I')

		self.Tree2.Branch('EF_e15_medium'   , self.EF_e15_medium   , 'EF_e15_medium/I')
		self.Tree2.Branch('EF_e20_medium'   , self.EF_e20_medium   , 'EF_e20_medium/I')
		self.Tree2.Branch('EF_e20_medium1'  , self.EF_e20_medium1  , 'EF_e20_medium1/I')
		self.Tree2.Branch('EF_2e12_medium'  , self.EF_2e12_medium  , 'EF_2e12_medium/I')

		self.Tree2.Branch('EF_2g20_loose'   , self.EF_2g20_loose   , 'EF_2g20_loose/I')

		self.Tree2.Branch('EF_mu10'         , self.EF_mu10         , 'EF_mu10/I')
		self.Tree2.Branch('EF_mu10_MG'      , self.EF_mu10_MG      , 'EF_mu10_MG/I')
		self.Tree2.Branch('EF_mu13_MG'      , self.EF_mu13_MG      , 'EF_mu13_MG/I')
		self.Tree2.Branch('EF_mu13_MG_tight', self.EF_mu13_MG_tight, 'EF_mu13_MG_tight/I')
		self.Tree2.Branch('EF_mu20_MG'      , self.EF_mu20_MG      , 'EF_mu20_MG/I')

	#####################################################################

	def treeCleaner(self):
		#########################
		# EVENT			#
		#########################

		self.RunNumber[0] = 0
		self.EventNumber[0] = 0
		self.lbn[0] = 0

		#########################
		# PRIMARY VERTICES	#
		#########################

		self.vxp_n[0] = 0

		self.vxp_x.clear()
		self.vxp_y.clear()
		self.vxp_z.clear()

		self.vxp_nTracks.clear()

		#########################
		# ELECTRONS		#
		#########################

		self.el_n[0] = 0

		self.el_E.clear()
		self.el_Et.clear()
		self.el_pt.clear()
		self.el_eta.clear()
		self.el_phi.clear()
		self.el_charge.clear()
		self.el_author.clear()
		self.el_isEM.clear()
		self.el_OQ.clear()

		self.el_cl_E.clear()
		self.el_cl_pt.clear()
		self.el_cl_eta.clear()
		self.el_cl_phi.clear()

		self.el_Es2.clear()
		self.el_etas2.clear()
		self.el_phis2.clear()

		self.el_loose.clear()
		self.el_medium.clear()
		self.el_tight.clear()

		self.el_reta.clear()
		self.el_weta2.clear()

		self.el_ptcone20.clear()
		self.el_ptcone30.clear()
		self.el_ptcone40.clear()

		self.el_Etcone20.clear()
		self.el_Etcone30.clear()
		self.el_Etcone40.clear()

		self.el_Etcone20_pt_corrected.clear()
		self.el_Etcone30_pt_corrected.clear()
		self.el_Etcone40_pt_corrected.clear()

		self.el_nBLHits.clear()
		self.el_nPixHits.clear()
		self.el_nSCTHits.clear()
		self.el_nTRTHits.clear()

		self.el_trackd0.clear()
		self.el_trackz0.clear()
		self.el_trackphi.clear()
		self.el_tracktheta.clear()
		self.el_trackqoverp.clear()

		self.el_trackd0pvunbiased.clear()
		self.el_trackz0pvunbiased.clear()
		self.el_tracksigd0pvunbiased.clear()
		self.el_tracksigz0pvunbiased.clear()

		self.el_truth_type.clear()
		self.el_truth_barbode.clear()
		self.el_truth_mothertype.clear()
		self.el_truth_motherbarcode.clear()

		self.el_type.clear()
		self.el_origin.clear()
		self.el_typebkg.clear()
		self.el_originbkg.clear()

		#########################
		# MUONS MUID		#
		#########################

		self.mu_muid_n[0] = 0

		self.mu_muid_E.clear()
		self.mu_muid_Et.clear()
		self.mu_muid_pt.clear()
		self.mu_muid_eta.clear()
		self.mu_muid_phi.clear()
		self.mu_muid_charge.clear()
		self.mu_muid_author.clear()

		self.mu_muid_loose.clear()
		self.mu_muid_medium.clear()
		self.mu_muid_tight.clear()

		self.mu_muid_ptcone20.clear()
		self.mu_muid_ptcone30.clear()
		self.mu_muid_ptcone40.clear()

		self.mu_muid_etcone20.clear()
		self.mu_muid_etcone30.clear()
		self.mu_muid_etcone40.clear()

		self.mu_muid_nBLHits.clear()
		self.mu_muid_nPixHits.clear()
		self.mu_muid_nSCTHits.clear()
		self.mu_muid_nTRTHits.clear()

		self.mu_muid_trackd0.clear()
		self.mu_muid_trackz0.clear()
		self.mu_muid_trackphi.clear()
		self.mu_muid_tracktheta.clear()
		self.mu_muid_trackqoverp.clear()

		self.mu_muid_trackd0pvunbiased.clear()
		self.mu_muid_trackz0pvunbiased.clear()
		self.mu_muid_tracksigd0pvunbiased.clear()
		self.mu_muid_tracksigz0pvunbiased.clear()

		self.mu_muid_truth_type.clear()
		self.mu_muid_truth_barbode.clear()
		self.mu_muid_truth_mothertype.clear()
		self.mu_muid_truth_motherbarcode.clear()

		#########################
		# MUONS STACO		#
		#########################

		self.mu_staco_n[0] = 0

		self.mu_staco_E.clear()
		self.mu_staco_Et.clear()
		self.mu_staco_pt.clear()
		self.mu_staco_eta.clear()
		self.mu_staco_phi.clear()
		self.mu_staco_charge.clear()
		self.mu_staco_author.clear()

		self.mu_staco_loose.clear()
		self.mu_staco_medium.clear()
		self.mu_staco_tight.clear()

		self.mu_staco_ptcone20.clear()
		self.mu_staco_ptcone30.clear()
		self.mu_staco_ptcone40.clear()

		self.mu_staco_etcone20.clear()
		self.mu_staco_etcone30.clear()
		self.mu_staco_etcone40.clear()

		self.mu_staco_nBLHits.clear()
		self.mu_staco_nPixHits.clear()
		self.mu_staco_nSCTHits.clear()
		self.mu_staco_nTRTHits.clear()

		self.mu_staco_trackd0.clear()
		self.mu_staco_trackz0.clear()
		self.mu_staco_trackphi.clear()
		self.mu_staco_tracktheta.clear()
		self.mu_staco_trackqoverp.clear()

		self.mu_staco_trackd0pvunbiased.clear()
		self.mu_staco_trackz0pvunbiased.clear()
		self.mu_staco_tracksigd0pvunbiased.clear()
		self.mu_staco_tracksigz0pvunbiased.clear()

		self.mu_staco_truth_type.clear()
		self.mu_staco_truth_barbode.clear()
		self.mu_staco_truth_mothertype.clear()
		self.mu_staco_truth_motherbarcode.clear()

		#########################
		# TRIGGERS		#
		#########################

		self.L1_EM14[0] = 0
		self.L1_MU10[0] = 0

		self.EF_e15_medium[0] = 0
		self.EF_e20_medium[0] = 0
		self.EF_e20_medium1[0] = 0
		self.EF_2e12_medium[0] = 0

		self.EF_2g20_loose[0] = 0

		self.EF_mu10[0] = 0
		self.EF_mu10_MG[0] = 0
		self.EF_mu13_MG[0] = 0
		self.EF_mu13_MG_tight[0] = 0
		self.EF_mu20_MG[0] = 0

	#####################################################################

	def finalize(self):
		return PyAthena.StatusCode.Success

	#####################################################################

	def execute(self):
		self.treeCleaner()

		#############################################################
		# EVENT							    #
		#############################################################

		L = self.StoreGateSvc.keys()

		if   'MyEvent' in L:
			event = self.StoreGateSvc['MyEvent']
		elif 'McEventInfo' in L:
			event = self.StoreGateSvc['McEventInfo']
		else:
			return PyAthena.StatusCode.Failure

		eventID = event.event_ID()

		#############################################################

		self.RunNumber[0] = eventID.run_number()
		self.EventNumber[0] = eventID.event_number()
		self.lbn[0] = eventID.lumi_block()

		#############################################################
		# PRIMARY VERTICES & PRIMARY TRACKS			    #
		#############################################################

		vertices = self.StoreGateSvc['VxPrimaryCandidate']

		#############################################################

		for vertex in vertices:
			self.vxp_x.push_back(vertex.recVertex().position().x())
			self.vxp_y.push_back(vertex.recVertex().position().y())
			self.vxp_z.push_back(vertex.recVertex().position().z())

			self.vxp_nTracks.push_back(len(vertex.vxTrackAtVertex()))

			##

			self.vxp_n[0] += 1

		#############################################################
		# ELECTRONS						    #
		#############################################################

		electrons = self.StoreGateSvc['ElectronAODCollection']

		#############################################################

		if len(vertices) > 0:
			for electron in electrons:

				E = electron.e()
				Et = electron.et()
				pt = electron.pt()
				eta = electron.eta()
				phi = electron.phi()
				charge = electron.charge()
				author = electron.author()
				isEM1 = ctypes.c_int32(electron.isem()).value
				isEM2 = ctypes.c_uint32(electron.isem()).value

				loose = False
				medium = False
				tight = False

				if electron.author() == 8:
					if electron.isElectron(PyAthena.egammaPID.frwdElectronLoose):
						loose = True
					if electron.isElectron(PyAthena.egammaPID.frwdElectronTight):
						tight = True
				else:
					if electron.isElectron(PyAthena.egammaPID.ElectronLoose):
						loose = True
					if electron.isElectron(PyAthena.egammaPID.ElectronMedium):
						medium = True
					if electron.isElectron(PyAthena.egammaPID.ElectronTight):
						tight = True

				##

				if electron.isgoodoq(PyAthena.egammaPID.BADCLUSELECTRON) == 0:
					OQ = 0x00000000000000000000000000000000
				else:
					OQ = PyAthena.egammaPID.BADCLUSELECTRON

				##

				is_ok = False

				for i in xrange(electron.nDetails()):

					detail = electron.detail(i)

					if detail and isinstance(detail, PyAthena.EMShower) and electron.detailName(i) == 'egDetailAOD':

						if detail.e277() != 0.0:
							reta = detail.e237() / detail.e277()
						else:
							reta =   -999999.0   /   +1.000000

						weta2 = detail.weta2()

						ptcone20 = detail.ptcone20()
						ptcone30 = detail.ptcone30()
						ptcone40 = detail.ptcone40()

						Etcone20 = detail.etcone20()
						Etcone30 = detail.etcone30()
						Etcone40 = detail.etcone40()

						Etcone20_pt_corrected = -999999.0
						Etcone30_pt_corrected = -999999.0
						Etcone40_pt_corrected = -999999.0

						is_ok = True

						break

				if is_ok == False:
					continue

				##

				cluster = electron.cluster()

				if cluster:
					cl_E = cluster.e()
					cl_pt = cluster.pt()
					cl_eta = cluster.eta()
					cl_phi = cluster.phi()

					cl_Es2 = cluster.energyBE(2)
					cl_etas2 = cluster.etaBE(2)
					cl_phis2 = cluster.phiBE(2)

				else:
					continue

				##

				track = electron.trackParticle()

				if track:
					summary = track.trackSummary()

					if summary:
						nBLHits = summary.get(ROOT.Trk.numberOfBLayerHits)
						nPixHits = summary.get(ROOT.Trk.numberOfPixelHits)
						nSCTHits = summary.get(ROOT.Trk.numberOfSCTHits)
						nTRTHits = summary.get(ROOT.Trk.numberOfTRTHits)
					else:
						continue

					perigee = track.measuredPerigee()

					if perigee:
						trackd0 = perigee.parameters()[0]
						trackz0 = perigee.parameters()[1]
						trackphi = perigee.parameters()[2]
						tracktheta = perigee.parameters()[3]
						trackqoverp = perigee.parameters()[4]
					else:
						continue

					r = self.TrackToVertexIPEstimator.estimate(track, vertices[0], True)

					trackd0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getD0(r)
					trackz0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getZ0(r)

					tracksigd0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getSigmaD0(r)
					tracksigz0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getSigmaZ0(r)

				else:
					continue

				##

				if isMC == False:
					truth_type = -999999
					truth_mothertype = -999999
					truth_barbode = -999999
					truth_motherbarcode = -999999

					epyt = -999999
					origin = -999999
					typebkg = -999999
					originbkg = -999999

				else:
					truth_type = -999999		# TODO #
					truth_mothertype = -999999	# TODO #
					truth_barbode = -999999		# TODO #
					truth_motherbarcode = -999999	# TODO #

					##

					epyt, origin = self.MCTruthClassifier.particleTruthClassifier(electron)

					if epyt > 0 \
					   and      \
					   epyt < 5:
						mc_electron = self.MCTruthClassifier.getGenPart()

						typebkg, originbkg = self.MCTruthClassifier.checkOrigOfBkgElec(mc_electron)
					else:
						typebkg, originbkg = -999999, -999999

				##

				self.el_E.push_back(E)
				self.el_Et.push_back(Et)
				self.el_pt.push_back(pt)
				self.el_eta.push_back(eta)
				self.el_phi.push_back(phi)
				self.el_charge.push_back(charge)
				self.el_author.push_back(author)
				self.el_isEM.push_back(isEM1)
				self.el_OQ.push_back(OQ)

				self.el_cl_E.push_back(cl_E)
				self.el_cl_pt.push_back(cl_pt)
				self.el_cl_eta.push_back(cl_eta)
				self.el_cl_phi.push_back(cl_phi)

				self.el_Es2.push_back(cl_Es2)
				self.el_etas2.push_back(cl_etas2)
				self.el_phis2.push_back(cl_phis2)

				self.el_loose.push_back(loose)
				self.el_medium.push_back(medium)
				self.el_tight.push_back(tight)

				self.el_reta.push_back(reta)
				self.el_weta2.push_back(weta2)

				self.el_ptcone20.push_back(ptcone20)
				self.el_ptcone30.push_back(ptcone30)
				self.el_ptcone40.push_back(ptcone40)

				self.el_Etcone20.push_back(Etcone20)
				self.el_Etcone30.push_back(Etcone30)
				self.el_Etcone40.push_back(Etcone40)

				self.el_Etcone20_pt_corrected.push_back(Etcone20_pt_corrected)
				self.el_Etcone30_pt_corrected.push_back(Etcone30_pt_corrected)
				self.el_Etcone40_pt_corrected.push_back(Etcone40_pt_corrected)

				self.el_nBLHits.push_back(nBLHits)
				self.el_nPixHits.push_back(nPixHits)
				self.el_nSCTHits.push_back(nSCTHits)
				self.el_nTRTHits.push_back(nTRTHits)

				self.el_trackd0.push_back(trackd0)
				self.el_trackz0.push_back(trackz0)
				self.el_trackphi.push_back(trackphi)
				self.el_tracktheta.push_back(tracktheta)
				self.el_trackqoverp.push_back(trackqoverp)

				self.el_trackd0pvunbiased.push_back(trackd0pvunbiased)
				self.el_trackz0pvunbiased.push_back(trackz0pvunbiased)

				self.el_tracksigd0pvunbiased.push_back(tracksigd0pvunbiased)
				self.el_tracksigz0pvunbiased.push_back(tracksigz0pvunbiased)

				self.el_truth_type.push_back(truth_type)
				self.el_truth_barbode.push_back(truth_barbode)
				self.el_truth_mothertype.push_back(truth_mothertype)
				self.el_truth_motherbarcode.push_back(truth_motherbarcode)

				self.el_type.push_back(epyt)
				self.el_origin.push_back(origin)
				self.el_typebkg.push_back(typebkg)
				self.el_originbkg.push_back(originbkg)

				##

				self.el_n[0] += 1

		#############################################################
		# MUONS MUID						    #
		#############################################################

		muons = self.StoreGateSvc['MuidMuonCollection']

		#############################################################

		if len(vertices) > 0:
			for muon in muons:

				E = muon.e()
				Et = muon.et()
				pt = muon.pt()
				eta = muon.eta()
				phi = muon.phi()
				charge = muon.charge()
				author = muon.author()

				##

				loose = muon.isLoose()
				medium = muon.isMedium()
				tight = muon.isTight()

				##

				ptcone20 = muon.parameter(PyAthena.MuonParameters.ptcone20)
				ptcone30 = muon.parameter(PyAthena.MuonParameters.ptcone30)
				ptcone40 = muon.parameter(PyAthena.MuonParameters.ptcone40)

				Etcone20 = muon.parameter(PyAthena.MuonParameters.etcone20)
				Etcone30 = muon.parameter(PyAthena.MuonParameters.etcone30)
				Etcone40 = muon.parameter(PyAthena.MuonParameters.etcone40)

				##

				track = muon.track()

				if track:
					summary = track.trackSummary()

					if summary:
						nBLHits = summary.get(ROOT.Trk.numberOfBLayerHits)
						nPixHits = summary.get(ROOT.Trk.numberOfPixelHits)
						nSCTHits = summary.get(ROOT.Trk.numberOfSCTHits)
						nTRTHits = summary.get(ROOT.Trk.numberOfTRTHits)
					else:
						continue

					perigee = track.measuredPerigee()

					if perigee:
						trackd0 = perigee.parameters()[0]
						trackz0 = perigee.parameters()[1]
						trackphi = perigee.parameters()[2]
						tracktheta = perigee.parameters()[3]
						trackqoverp = perigee.parameters()[4]
					else:
						continue

					r = self.TrackToVertexIPEstimator.estimate(track, vertices[0], True)

					trackd0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getD0(r)
					trackz0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getZ0(r)

					tracksigd0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getSigmaD0(r)
					tracksigz0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getSigmaZ0(r)

				else:
					continue

				##

				if isMC == False:
					truth_type = -999999
					truth_mothertype = -999999
					truth_barbode = -999999
					truth_motherbarcode = -999999

				else:
					truth_type = -999999		# TODO #
					truth_mothertype = -999999	# TODO #
					truth_barbode = -999999		# TODO #
					truth_motherbarcode = -999999	# TODO #

				##

				self.mu_muid_E.push_back(E)
				self.mu_muid_Et.push_back(Et)
				self.mu_muid_pt.push_back(pt)
				self.mu_muid_eta.push_back(eta)
				self.mu_muid_phi.push_back(phi)
				self.mu_muid_charge.push_back(charge)
				self.mu_muid_author.push_back(author)

				self.mu_muid_loose.push_back(loose)
				self.mu_muid_medium.push_back(medium)
				self.mu_muid_tight.push_back(tight)

				self.mu_muid_ptcone20.push_back(ptcone20)
				self.mu_muid_ptcone30.push_back(ptcone30)
				self.mu_muid_ptcone40.push_back(ptcone40)

				self.mu_muid_etcone20.push_back(Etcone20)
				self.mu_muid_etcone30.push_back(Etcone30)
				self.mu_muid_etcone40.push_back(Etcone40)

				self.mu_muid_nBLHits.push_back(nBLHits)
				self.mu_muid_nPixHits.push_back(nPixHits)
				self.mu_muid_nSCTHits.push_back(nSCTHits)
				self.mu_muid_nTRTHits.push_back(nTRTHits)

				self.mu_muid_trackd0.push_back(trackd0)
				self.mu_muid_trackz0.push_back(trackz0)
				self.mu_muid_trackphi.push_back(trackphi)
				self.mu_muid_tracktheta.push_back(tracktheta)
				self.mu_muid_trackqoverp.push_back(trackqoverp)

				self.mu_muid_trackd0pvunbiased.push_back(trackd0pvunbiased)
				self.mu_muid_trackz0pvunbiased.push_back(trackz0pvunbiased)
				self.mu_muid_tracksigd0pvunbiased.push_back(tracksigd0pvunbiased)
				self.mu_muid_tracksigz0pvunbiased.push_back(tracksigz0pvunbiased)

				self.mu_muid_truth_type.push_back(truth_type)
				self.mu_muid_truth_barbode.push_back(truth_barbode)
				self.mu_muid_truth_mothertype.push_back(truth_mothertype)
				self.mu_muid_truth_motherbarcode.push_back(truth_motherbarcode)

				##

				self.mu_muid_n[0] += 1

		#############################################################
		# MUONS STACO						    #
		#############################################################

		muons = self.StoreGateSvc['StacoMuonCollection']

		#############################################################

		if len(vertices) > 0:
			for muon in muons:

				E = muon.e()
				Et = muon.et()
				pt = muon.pt()
				eta = muon.eta()
				phi = muon.phi()
				charge = muon.charge()
				author = muon.author()

				##

				loose = muon.isLoose()
				medium = muon.isMedium()
				tight = muon.isTight()

				##

				ptcone20 = muon.parameter(PyAthena.MuonParameters.ptcone20)
				ptcone30 = muon.parameter(PyAthena.MuonParameters.ptcone30)
				ptcone40 = muon.parameter(PyAthena.MuonParameters.ptcone40)

				Etcone20 = muon.parameter(PyAthena.MuonParameters.etcone20)
				Etcone30 = muon.parameter(PyAthena.MuonParameters.etcone30)
				Etcone40 = muon.parameter(PyAthena.MuonParameters.etcone40)

				##

				track = muon.track()

				if track:
					summary = track.trackSummary()

					if summary:
						nBLHits = summary.get(ROOT.Trk.numberOfBLayerHits)
						nPixHits = summary.get(ROOT.Trk.numberOfPixelHits)
						nSCTHits = summary.get(ROOT.Trk.numberOfSCTHits)
						nTRTHits = summary.get(ROOT.Trk.numberOfTRTHits)
					else:
						continue

					perigee = track.measuredPerigee()

					if perigee:
						trackd0 = perigee.parameters()[0]
						trackz0 = perigee.parameters()[1]
						trackphi = perigee.parameters()[2]
						tracktheta = perigee.parameters()[3]
						trackqoverp = perigee.parameters()[4]
					else:
						continue

					r = self.TrackToVertexIPEstimator.estimate(track, vertices[0], True)

					trackd0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getD0(r)
					trackz0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getZ0(r)

					tracksigd0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getSigmaD0(r)
					tracksigz0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getSigmaZ0(r)

				else:
					continue

				##

				if isMC == False:
					truth_type = -999999
					truth_mothertype = -999999
					truth_barbode = -999999
					truth_motherbarcode = -999999

				else:
					truth_type = -999999		# TODO #
					truth_mothertype = -999999	# TODO #
					truth_barbode = -999999		# TODO #
					truth_motherbarcode = -999999	# TODO #

				##

				self.mu_staco_E.push_back(E)
				self.mu_staco_Et.push_back(Et)
				self.mu_staco_pt.push_back(pt)
				self.mu_staco_eta.push_back(eta)
				self.mu_staco_phi.push_back(phi)
				self.mu_staco_charge.push_back(charge)
				self.mu_staco_author.push_back(author)

				self.mu_staco_loose.push_back(loose)
				self.mu_staco_medium.push_back(medium)
				self.mu_staco_tight.push_back(tight)

				self.mu_staco_ptcone20.push_back(ptcone20)
				self.mu_staco_ptcone30.push_back(ptcone30)
				self.mu_staco_ptcone40.push_back(ptcone40)

				self.mu_staco_etcone20.push_back(Etcone20)
				self.mu_staco_etcone30.push_back(Etcone30)
				self.mu_staco_etcone40.push_back(Etcone40)

				self.mu_staco_nBLHits.push_back(nBLHits)
				self.mu_staco_nPixHits.push_back(nPixHits)
				self.mu_staco_nSCTHits.push_back(nSCTHits)
				self.mu_staco_nTRTHits.push_back(nTRTHits)

				self.mu_staco_trackd0.push_back(trackd0)
				self.mu_staco_trackz0.push_back(trackz0)
				self.mu_staco_trackphi.push_back(trackphi)
				self.mu_staco_tracktheta.push_back(tracktheta)
				self.mu_staco_trackqoverp.push_back(trackqoverp)

				self.mu_staco_trackd0pvunbiased.push_back(trackd0pvunbiased)
				self.mu_staco_trackz0pvunbiased.push_back(trackz0pvunbiased)
				self.mu_staco_tracksigd0pvunbiased.push_back(tracksigd0pvunbiased)
				self.mu_staco_tracksigz0pvunbiased.push_back(tracksigz0pvunbiased)

				self.mu_staco_truth_type.push_back(truth_type)
				self.mu_staco_truth_barbode.push_back(truth_barbode)
				self.mu_staco_truth_mothertype.push_back(truth_mothertype)
				self.mu_staco_truth_motherbarcode.push_back(truth_motherbarcode)

				##

				self.mu_staco_n[0] += 1

		#############################################################
		# TRIGGERS						    #
		#############################################################

		ctp = self.StoreGateSvc['CTP_Decision']

		ctp_items = ctp.getItems()

		#############################################################

		self.L1_EM14[0] = 0
		self.L1_MU10[0] = 0

		for ctp_item in ctp_items:

			if ctp_item == 'L1_EM14':
				self.L1_EM14[0] = 1

			if ctp_item == 'L1_MU10':
				self.L1_MU10[0] = 1

		##

		if self.TrigDecisionTool.getChainGroup('EF_e15_medium').isPassed():
			self.EF_e15_medium[0] = 1
		else:
			self.EF_e15_medium[0] = 0

		if self.TrigDecisionTool.getChainGroup('EF_e20_medium').isPassed():
			self.EF_e20_medium[0] = 1
		else:
			self.EF_e20_medium[0] = 0

		if self.TrigDecisionTool.getChainGroup('EF_e20_medium1').isPassed():
			self.EF_e20_medium1[0] = 1
		else:
			self.EF_e20_medium1[0] = 0

		if self.TrigDecisionTool.getChainGroup('EF_2e12_medium').isPassed():
			self.EF_2e12_medium[0] = 1
		else:
			self.EF_2e12_medium[0] = 0

		##

		if self.TrigDecisionTool.getChainGroup('EF_2e12_medium').isPassed():
			self.EF_2g20_loose[0] = 1
		else:
			self.EF_2g20_loose[0] = 0

		##

		if self.TrigDecisionTool.getChainGroup('EF_mu10').isPassed():
			self.EF_mu10[0] = 1
		else:
			self.EF_mu10[0] = 0

		if self.TrigDecisionTool.getChainGroup('EF_mu10_MG').isPassed():
			self.EF_mu10_MG[0] = 1
		else:
			self.EF_mu10_MG[0] = 0

		if self.TrigDecisionTool.getChainGroup('EF_mu13_MG').isPassed():
			self.EF_mu13_MG[0] = 1
		else:
			self.EF_mu13_MG[0] = 0

		if self.TrigDecisionTool.getChainGroup('EF_mu13_MG_tight').isPassed():
			self.EF_mu13_MG_tight[0] = 1
		else:
			self.EF_mu13_MG_tight[0] = 0

		if self.TrigDecisionTool.getChainGroup('EF_mu20_MG').isPassed():
			self.EF_mu20_MG[0] = 1
		else:
			self.EF_mu20_MG[0] = 0

		#############################################################
		#############################################################

		self.Tree1.Fill()
		if isEGamma:
			self.Tree2.Fill()

		return PyAthena.StatusCode.Success

#############################################################################

from AthenaCommon.AlgSequence import AlgSequence

job = AlgSequence()
job += uD3PD('uD3PD')

#############################################################################

