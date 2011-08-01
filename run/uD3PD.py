#############################################################################
# USER OPTIONS								    #
#############################################################################

#AtlGeo = 'ATLAS-GEO-16-00-00'
#CondDB = 'OFLCOND-SDR-BS7T-04-13'

AtlGeo = 'ATLAS-GEO-16-00-01'
CondDB = 'COMCOND-BLKPST-004-07'

#############################################################################

isMC = False

isEGamma = False

InputFormat = 'AOD'

#############################################################################

InputFiles = [
#	'../../AOD.280342._000152.pool.root'
	'/tmp/jodier/DAOD_HSG2.384735._000027.pool.root'
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

theTrigDecisionTool = Trig__TrigDecisionTool(name = 'TrigDecisionTool')
theTrigDecisionTool.OutputLevel = ERROR

ToolSvc += theTrigDecisionTool

#################################
# AtlasExtrapolator		#
#################################

from TrkExTools.AtlasExtrapolator import AtlasExtrapolator

theAtlasExtrapolator = AtlasExtrapolator(name = 'AtlasExtrapolator')
theAtlasExtrapolator.OutputLevel = ERROR

ToolSvc += theAtlasExtrapolator

#################################
# ExtrapolateToCaloTool		#
#################################

from TrackToCalo.TrackToCaloConf import ExtrapolateToCaloTool

theExtrapolateToCaloTool = ExtrapolateToCaloTool(name = 'ExtrapolateToCaloTool', Extrapolator = theAtlasExtrapolator)
theExtrapolateToCaloTool.OutputLevel = ERROR

ToolSvc += theExtrapolateToCaloTool

#################################
# MCTruthClassifier		#
#################################

from MCTruthClassifier.MCTruthClassifierConf import MCTruthClassifier

theMCTruthClassifier = MCTruthClassifier(name = 'MCTruthClassifier', ExtrapolateToCaloTool = theExtrapolateToCaloTool)
theMCTruthClassifier.OutputLevel = ERROR

ToolSvc += theMCTruthClassifier

#################################
# TrackToVertexIPEstimator	#
#################################

from TrkVertexFitterUtils.TrkVertexFitterUtilsConf import Trk__TrackToVertexIPEstimator

theTrackToVertexIPEstimator = Trk__TrackToVertexIPEstimator(name = 'TrackToVertexIPEstimator', Extrapolator = theAtlasExtrapolator)
theTrackToVertexIPEstimator.OutputLevel = ERROR

ToolSvc += theTrackToVertexIPEstimator

#############################################################################
# ALGORITHM								    #
#############################################################################

import ROOT
import PyCintex

#############################################################################

import AthenaPython.PyAthena as PyAthena

#############################################################################

import math
import array
import ctypes

#############################################################################

def __dR2(eta1, eta2, phi1, phi2):

	dEta = eta1 - eta2
	dPhi = phi1 - phi2

	while dPhi < -math.pi:
		dPhi += 2.0 * math.pi

	while dPhi >= +math.pi:
		dPhi -= 2.0 * math.pi

	return dEta * dEta + dPhi * dPhi

#############################################################################

def particleMatching(theEta, thePhi, container, radius = 0.15):

	goodI = 0
	goodDR = 999999.0

	for i in xrange(len(container)):

		#############################################################

		T = type(container[i]).__name__

		#############################################################

		if   T == 'TrigMuonEFInfo':
			if container[i].hasExtrapolatedTrack():
				eta = container[i].ExtrapolatedTrack().eta()
				phi = container[i].ExtrapolatedTrack().phi()
			else:
				eta = -999999.0
				phi = -999999.0
		elif T == 'Trig::Feature<TrigRoiDescriptor>':
			eta = container[i].cptr().eta0()
			phi = container[i].cptr().phi0()
		else:
			eta = container[i].eta()
			phi = container[i].phi()

		#############################################################

		dR = math.sqrt(__dR2(theEta, eta, thePhi, phi))

		if goodDR > dR:
			goodI = i
			goodDR = dR

		#############################################################

	if goodDR < radius:
		return [goodI, goodDR]
	else:
		return [  -1 ,   -1  ]

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

		if not self.TrigDecisionTool:
			return PyAthena.StatusCode.Failure

		##

		self.MCTruthClassifier = PyAthena.py_tool('MCTruthClassifier', iface = 'IMCTruthClassifier')

		if not self.MCTruthClassifier:
			return PyAthena.StatusCode.Failure

		##

		self.TrackToVertexIPEstimator = PyAthena.py_tool('Trk::TrackToVertexIPEstimator/TrackToVertexIPEstimator', iface = 'Trk::ITrackToVertexIPEstimator')

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

		self.pixelError = array.array('I', [0])
		self.sctError = array.array('I', [0])
		self.trtError = array.array('I', [0])
		self.larError = array.array('I', [0])
		self.muonError = array.array('I', [0])

		self.mcevt_weight = ROOT.std.vector(float)()

		##

		self.Tree1.Branch('RunNumber', self.RunNumber, 'RunNumber/i')
		self.Tree1.Branch('EventNumber', self.EventNumber, 'EventNumber/i')
		self.Tree1.Branch('lbn', self.lbn, 'lbn/i')

		self.Tree1.Branch('pixelError', self.pixelError, 'pixelError/i')
		self.Tree1.Branch('sctError', self.sctError, 'sctError/i')
		self.Tree1.Branch('trtError', self.trtError, 'trtError/i')
		self.Tree1.Branch('larError', self.larError, 'larError/i')
		self.Tree1.Branch('muonError', self.muonError, 'muonError/i')

		self.Tree1.Branch('mcevt_weight', self.mcevt_weight)

		if isEGamma:
			self.Tree2.Branch('RunNumber', self.RunNumber, 'RunNumber/i')
			self.Tree2.Branch('EventNumber', self.EventNumber, 'EventNumber/i')
			self.Tree2.Branch('lbn', self.lbn, 'lbn/i')

			self.Tree2.Branch('pixelError', self.pixelError, 'pixelError/i')
			self.Tree2.Branch('sctError', self.sctError, 'sctError/i')
			self.Tree2.Branch('trtError', self.trtError, 'trtError/i')
			self.Tree2.Branch('larError', self.larError, 'larError/i')
			self.Tree2.Branch('muonError', self.muonError, 'muonError/i')

			self.Tree2.Branch('mcevt_weight', self.mcevt_weight)

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

		self.el_m = ROOT.std.vector(float)()
		self.el_E = ROOT.std.vector(float)()
		self.el_Et = ROOT.std.vector(float)()
		self.el_pt = ROOT.std.vector(float)()
		self.el_eta = ROOT.std.vector(float)()
		self.el_phi = ROOT.std.vector(float)()
		self.el_charge = ROOT.std.vector(float)()
		self.el_author = ROOT.std.vector(int)()
		self.el_isEM = ROOT.std.vector("unsigned int")()
		self.el_OQ = ROOT.std.vector("unsigned int")()

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

		self.el_E237 = ROOT.std.vector(float)()
		self.el_E277 = ROOT.std.vector(float)()
		self.el_emaxs1 = ROOT.std.vector(float)()
		self.el_Emax2 = ROOT.std.vector(float)()
		self.el_Ethad = ROOT.std.vector(float)()
		self.el_Ethad1 = ROOT.std.vector(float)()
		self.el_f1 = ROOT.std.vector(float)()
		self.el_weta2 = ROOT.std.vector(float)()
		self.el_wstot = ROOT.std.vector(float)()
		self.el_etap = ROOT.std.vector(float)()

		self.el_nBLHits = ROOT.std.vector(int)()
		self.el_nPixHits = ROOT.std.vector(int)()
		self.el_nSCTHits = ROOT.std.vector(int)()
		self.el_nTRTHits = ROOT.std.vector(int)()
		self.el_nTRTHighTHits = ROOT.std.vector(int)()

		self.el_nBLOutliers = ROOT.std.vector(int)()
		self.el_nPixOutliers = ROOT.std.vector(int)()
		self.el_nSCTOutliers = ROOT.std.vector(int)()
		self.el_nTRTOutliers = ROOT.std.vector(int)()
		self.el_nTRTHighTOutliers = ROOT.std.vector(int)()

		self.el_nPixHoles = ROOT.std.vector(int)()
		self.el_nSCTHoles = ROOT.std.vector(int)()
		self.el_nTRTHoles = ROOT.std.vector(int)()

		self.el_expectBLayerHit = ROOT.std.vector(int)()

		self.el_deltaeta1 = ROOT.std.vector(float)()
		self.el_deltaeta2 = ROOT.std.vector(float)()

		self.el_trackd0 = ROOT.std.vector(float)()
		self.el_trackz0 = ROOT.std.vector(float)()
		self.el_trackpt = ROOT.std.vector(float)()
		self.el_tracketa = ROOT.std.vector(float)()
		self.el_trackphi = ROOT.std.vector(float)()
		self.el_tracktheta = ROOT.std.vector(float)()
		self.el_trackqoverp = ROOT.std.vector(float)()

		self.el_trackd0pvunbiased = ROOT.std.vector(float)()
		self.el_trackz0pvunbiased = ROOT.std.vector(float)()
		self.el_tracksigd0pvunbiased = ROOT.std.vector(float)()
		self.el_tracksigz0pvunbiased = ROOT.std.vector(float)()

		self.el_ptcone20 = ROOT.std.vector(float)()
		self.el_ptcone30 = ROOT.std.vector(float)()
		self.el_ptcone40 = ROOT.std.vector(float)()

		self.el_Etcone20 = ROOT.std.vector(float)()
		self.el_Etcone30 = ROOT.std.vector(float)()
		self.el_Etcone40 = ROOT.std.vector(float)()

		self.el_truth_type = ROOT.std.vector(int)()
		self.el_truth_mothertype = ROOT.std.vector(int)()
		self.el_truth_barcode = ROOT.std.vector(int)()
		self.el_truth_motherbarcode = ROOT.std.vector(int)()

		self.el_type = ROOT.std.vector(int)()
		self.el_origin = ROOT.std.vector(int)()
		self.el_typebkg = ROOT.std.vector(int)()
		self.el_originbkg = ROOT.std.vector(int)()

		self.el_EF_dr = ROOT.std.vector(float)()
		self.el_EF_index = ROOT.std.vector(int)()

		##

		self.Tree1.Branch('el_n', self.el_n, 'el_n/I')

		self.Tree1.Branch('el_m', self.el_m)
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

		self.Tree1.Branch('el_E237', self.el_E237)
		self.Tree1.Branch('el_E277', self.el_E277)
		self.Tree1.Branch('el_emaxs1', self.el_emaxs1)
		self.Tree1.Branch('el_Emax2', self.el_Emax2)
		self.Tree1.Branch('el_Ethad', self.el_Ethad)
		self.Tree1.Branch('el_Ethad1', self.el_Ethad1)
		self.Tree1.Branch('el_f1', self.el_f1)
		self.Tree1.Branch('el_weta2', self.el_weta2)
		self.Tree1.Branch('el_wstot', self.el_wstot)
		self.Tree1.Branch('el_etap', self.el_etap)

		self.Tree1.Branch('el_nBLHits', self.el_nBLHits)
		self.Tree1.Branch('el_nPixHits', self.el_nPixHits)
		self.Tree1.Branch('el_nSCTHits', self.el_nSCTHits)
		self.Tree1.Branch('el_nTRTHits', self.el_nTRTHits)
		self.Tree1.Branch('el_nTRTHighTHits', self.el_nTRTHighTHits)

		self.Tree1.Branch('el_nBLOutliers', self.el_nBLOutliers)
		self.Tree1.Branch('el_nPixOutliers', self.el_nPixOutliers)
		self.Tree1.Branch('el_nSCTOutliers', self.el_nSCTOutliers)
		self.Tree1.Branch('el_nTRTOutliers', self.el_nTRTOutliers)
		self.Tree1.Branch('el_nTRTHighTOutliers', self.el_nTRTHighTOutliers)

		self.Tree1.Branch('el_nPixHoles', self.el_nPixHoles)
		self.Tree1.Branch('el_nSCTHoles', self.el_nSCTHoles)
		self.Tree1.Branch('el_nTRTHoles', self.el_nTRTHoles)

		self.Tree1.Branch('el_expectBLayerHit', self.el_expectBLayerHit)

		self.Tree1.Branch('el_deltaeta1', self.el_deltaeta1)
		self.Tree1.Branch('el_deltaeta2', self.el_deltaeta2)

		self.Tree1.Branch('el_trackd0', self.el_trackd0)
		self.Tree1.Branch('el_trackz0', self.el_trackz0)
		self.Tree1.Branch('el_trackpt', self.el_trackpt)
		self.Tree1.Branch('el_tracketa', self.el_tracketa)
		self.Tree1.Branch('el_trackphi', self.el_trackphi)
		self.Tree1.Branch('el_tracktheta', self.el_tracktheta)
		self.Tree1.Branch('el_trackqoverp', self.el_trackqoverp)

		self.Tree1.Branch('el_trackd0pvunbiased', self.el_trackd0pvunbiased)
		self.Tree1.Branch('el_trackz0pvunbiased', self.el_trackz0pvunbiased)
		self.Tree1.Branch('el_tracksigd0pvunbiased', self.el_tracksigd0pvunbiased)
		self.Tree1.Branch('el_tracksigz0pvunbiased', self.el_tracksigz0pvunbiased)

		self.Tree1.Branch('el_ptcone20', self.el_ptcone20)
		self.Tree1.Branch('el_ptcone30', self.el_ptcone30)
		self.Tree1.Branch('el_ptcone40', self.el_ptcone40)

		self.Tree1.Branch('el_Etcone20', self.el_Etcone20)
		self.Tree1.Branch('el_Etcone30', self.el_Etcone30)
		self.Tree1.Branch('el_Etcone40', self.el_Etcone40)

		self.Tree1.Branch('el_truth_type', self.el_truth_type)
		self.Tree1.Branch('el_truth_mothertype', self.el_truth_mothertype)
		self.Tree1.Branch('el_truth_barcode', self.el_truth_barcode)
		self.Tree1.Branch('el_truth_motherbarcode', self.el_truth_motherbarcode)

		self.Tree1.Branch('el_type', self.el_type)
		self.Tree1.Branch('el_origin', self.el_origin)
		self.Tree1.Branch('el_typebkg', self.el_typebkg)
		self.Tree1.Branch('el_originbkg', self.el_originbkg)

		self.Tree1.Branch('el_EF_dr', self.el_EF_dr)
		self.Tree1.Branch('el_EF_index', self.el_EF_index)

		#########################
		# MUONS MUID		#
		#########################

		self.mu_muid_n = array.array('i', [0])

		self.mu_muid_m = ROOT.std.vector(float)()
		self.mu_muid_E = ROOT.std.vector(float)()
		self.mu_muid_pt = ROOT.std.vector(float)()
		self.mu_muid_eta = ROOT.std.vector(float)()
		self.mu_muid_phi = ROOT.std.vector(float)()
		self.mu_muid_charge = ROOT.std.vector(float)()
		self.mu_muid_author = ROOT.std.vector(int)()

		self.mu_muid_loose = ROOT.std.vector(int)()
		self.mu_muid_medium = ROOT.std.vector(int)()
		self.mu_muid_tight = ROOT.std.vector(int)()

		self.mu_muid_nBLHits = ROOT.std.vector(int)()
		self.mu_muid_nPixHits = ROOT.std.vector(int)()
		self.mu_muid_nSCTHits = ROOT.std.vector(int)()
		self.mu_muid_nTRTHits = ROOT.std.vector(int)()
		self.mu_muid_nTRTHighTHits = ROOT.std.vector(int)()

		self.mu_muid_nBLOutliers = ROOT.std.vector(int)()
		self.mu_muid_nPixOutliers = ROOT.std.vector(int)()
		self.mu_muid_nSCTOutliers = ROOT.std.vector(int)()
		self.mu_muid_nTRTOutliers = ROOT.std.vector(int)()
		self.mu_muid_nTRTHighTOutliers = ROOT.std.vector(int)()

		self.mu_muid_nPixHoles = ROOT.std.vector(int)()
		self.mu_muid_nSCTHoles = ROOT.std.vector(int)()
		self.mu_muid_nTRTHoles = ROOT.std.vector(int)()

		self.mu_muid_expectBLayerHit = ROOT.std.vector(int)()

		self.mu_muid_nPixDeadSensors = ROOT.std.vector(int)()
		self.mu_muid_nSCTDeadSensors = ROOT.std.vector(int)()

		self.mu_muid_id_d0 = ROOT.std.vector(float)()
		self.mu_muid_id_z0 = ROOT.std.vector(float)()
		self.mu_muid_id_phi = ROOT.std.vector(float)()
		self.mu_muid_id_theta = ROOT.std.vector(float)()
		self.mu_muid_id_qoverp = ROOT.std.vector(float)()

		self.mu_muid_id_theta_exPV = ROOT.std.vector(float)()
		self.mu_muid_id_qoverp_exPV = ROOT.std.vector(float)()

		self.mu_muid_trackd0pvunbiased = ROOT.std.vector(float)()
		self.mu_muid_trackz0pvunbiased = ROOT.std.vector(float)()
		self.mu_muid_tracksigd0pvunbiased = ROOT.std.vector(float)()
		self.mu_muid_tracksigz0pvunbiased = ROOT.std.vector(float)()

		self.mu_muid_ptcone20 = ROOT.std.vector(float)()
		self.mu_muid_ptcone30 = ROOT.std.vector(float)()
		self.mu_muid_ptcone40 = ROOT.std.vector(float)()

		self.mu_muid_etcone20 = ROOT.std.vector(float)()
		self.mu_muid_etcone30 = ROOT.std.vector(float)()
		self.mu_muid_etcone40 = ROOT.std.vector(float)()

		self.mu_muid_truth_type = ROOT.std.vector(int)()
		self.mu_muid_truth_mothertype = ROOT.std.vector(int)()
		self.mu_muid_truth_barcode = ROOT.std.vector(int)()
		self.mu_muid_truth_motherbarcode = ROOT.std.vector(int)()

		self.mu_muid_EFCB_dr = ROOT.std.vector(float)()
		self.mu_muid_EFCB_index = ROOT.std.vector(int)()

		##

		self.Tree1.Branch('mu_muid_n', self.mu_muid_n, 'mu_muid_n/I')

		self.Tree1.Branch('mu_muid_m', self.mu_muid_m)
		self.Tree1.Branch('mu_muid_E', self.mu_muid_E)
		self.Tree1.Branch('mu_muid_pt', self.mu_muid_pt)
		self.Tree1.Branch('mu_muid_eta', self.mu_muid_eta)
		self.Tree1.Branch('mu_muid_phi', self.mu_muid_phi)
		self.Tree1.Branch('mu_muid_charge', self.mu_muid_charge)
		self.Tree1.Branch('mu_muid_author', self.mu_muid_author)

		self.Tree1.Branch('mu_muid_loose', self.mu_muid_loose)
		self.Tree1.Branch('mu_muid_medium', self.mu_muid_medium)
		self.Tree1.Branch('mu_muid_tight', self.mu_muid_tight)

		self.Tree1.Branch('mu_muid_nBLHits', self.mu_muid_nBLHits)
		self.Tree1.Branch('mu_muid_nPixHits', self.mu_muid_nPixHits)
		self.Tree1.Branch('mu_muid_nSCTHits', self.mu_muid_nSCTHits)
		self.Tree1.Branch('mu_muid_nTRTHits', self.mu_muid_nTRTHits)
		self.Tree1.Branch('mu_muid_nTRTHighTHits', self.mu_muid_nTRTHighTHits)

		self.Tree1.Branch('mu_muid_nBLOutliers', self.mu_muid_nBLOutliers)
		self.Tree1.Branch('mu_muid_nPixOutliers', self.mu_muid_nPixOutliers)
		self.Tree1.Branch('mu_muid_nSCTOutliers', self.mu_muid_nSCTOutliers)
		self.Tree1.Branch('mu_muid_nTRTOutliers', self.mu_muid_nTRTOutliers)
		self.Tree1.Branch('mu_muid_nTRTHighTOutliers', self.mu_muid_nTRTHighTOutliers)

		self.Tree1.Branch('mu_muid_nPixHoles', self.mu_muid_nPixHoles)
		self.Tree1.Branch('mu_muid_nSCTHoles', self.mu_muid_nSCTHoles)
		self.Tree1.Branch('mu_muid_nTRTHoles', self.mu_muid_nTRTHoles)

		self.Tree1.Branch('mu_muid_expectBLayerHit', self.mu_muid_expectBLayerHit)

		self.Tree1.Branch('mu_muid_nPixDeadSensors', self.mu_muid_nPixDeadSensors)
		self.Tree1.Branch('mu_muid_nSCTDeadSensors', self.mu_muid_nSCTDeadSensors)

		self.Tree1.Branch('mu_muid_id_d0', self.mu_muid_id_d0)
		self.Tree1.Branch('mu_muid_id_z0', self.mu_muid_id_z0)
		self.Tree1.Branch('mu_muid_id_phi', self.mu_muid_id_phi)
		self.Tree1.Branch('mu_muid_id_theta', self.mu_muid_id_theta)
		self.Tree1.Branch('mu_muid_id_qoverp', self.mu_muid_id_qoverp)

		self.Tree1.Branch('mu_muid_id_theta_exPV', self.mu_muid_id_theta_exPV)
		self.Tree1.Branch('mu_muid_id_qoverp_exPV', self.mu_muid_id_qoverp_exPV)

		self.Tree1.Branch('mu_muid_trackd0pvunbiased', self.mu_muid_trackd0pvunbiased)
		self.Tree1.Branch('mu_muid_trackz0pvunbiased', self.mu_muid_trackz0pvunbiased)
		self.Tree1.Branch('mu_muid_tracksigd0pvunbiased', self.mu_muid_tracksigd0pvunbiased)
		self.Tree1.Branch('mu_muid_tracksigz0pvunbiased', self.mu_muid_tracksigz0pvunbiased)

		self.Tree1.Branch('mu_muid_ptcone20', self.mu_muid_ptcone20)
		self.Tree1.Branch('mu_muid_ptcone30', self.mu_muid_ptcone30)
		self.Tree1.Branch('mu_muid_ptcone40', self.mu_muid_ptcone40)

		self.Tree1.Branch('mu_muid_etcone20', self.mu_muid_etcone20)
		self.Tree1.Branch('mu_muid_etcone30', self.mu_muid_etcone30)
		self.Tree1.Branch('mu_muid_etcone40', self.mu_muid_etcone40)

		self.Tree1.Branch('mu_muid_truth_type', self.mu_muid_truth_type)
		self.Tree1.Branch('mu_muid_truth_mothertype', self.mu_muid_truth_mothertype)
		self.Tree1.Branch('mu_muid_truth_barcode', self.mu_muid_truth_barcode)
		self.Tree1.Branch('mu_muid_truth_motherbarcode', self.mu_muid_truth_motherbarcode)

		self.Tree1.Branch('mu_muid_EFCB_dr', self.mu_muid_EFCB_dr)
		self.Tree1.Branch('mu_muid_EFCB_index', self.mu_muid_EFCB_index)

		#########################
		# MUONS STACO		#
		#########################

		self.mu_staco_n = array.array('i', [0])

		self.mu_staco_m = ROOT.std.vector(float)()
		self.mu_staco_E = ROOT.std.vector(float)()
		self.mu_staco_pt = ROOT.std.vector(float)()
		self.mu_staco_eta = ROOT.std.vector(float)()
		self.mu_staco_phi = ROOT.std.vector(float)()
		self.mu_staco_charge = ROOT.std.vector(float)()
		self.mu_staco_author = ROOT.std.vector(int)()

		self.mu_staco_loose = ROOT.std.vector(int)()
		self.mu_staco_medium = ROOT.std.vector(int)()
		self.mu_staco_tight = ROOT.std.vector(int)()

		self.mu_staco_nBLHits = ROOT.std.vector(int)()
		self.mu_staco_nPixHits = ROOT.std.vector(int)()
		self.mu_staco_nSCTHits = ROOT.std.vector(int)()
		self.mu_staco_nTRTHits = ROOT.std.vector(int)()
		self.mu_staco_nTRTHighTHits = ROOT.std.vector(int)()

		self.mu_staco_nBLHits = ROOT.std.vector(int)()
		self.mu_staco_nPixHits = ROOT.std.vector(int)()
		self.mu_staco_nSCTHits = ROOT.std.vector(int)()
		self.mu_staco_nTRTHits = ROOT.std.vector(int)()
		self.mu_staco_nTRTHighTHits = ROOT.std.vector(int)()

		self.mu_staco_nBLOutliers = ROOT.std.vector(int)()
		self.mu_staco_nPixOutliers = ROOT.std.vector(int)()
		self.mu_staco_nSCTOutliers = ROOT.std.vector(int)()
		self.mu_staco_nTRTOutliers = ROOT.std.vector(int)()
		self.mu_staco_nTRTHighTOutliers = ROOT.std.vector(int)()

		self.mu_staco_nPixHoles = ROOT.std.vector(int)()
		self.mu_staco_nSCTHoles = ROOT.std.vector(int)()
		self.mu_staco_nTRTHoles = ROOT.std.vector(int)()

		self.mu_staco_expectBLayerHit = ROOT.std.vector(int)()

		self.mu_staco_nPixDeadSensors = ROOT.std.vector(int)()
		self.mu_staco_nSCTDeadSensors = ROOT.std.vector(int)()

		self.mu_staco_id_d0 = ROOT.std.vector(float)()
		self.mu_staco_id_z0 = ROOT.std.vector(float)()
		self.mu_staco_id_phi = ROOT.std.vector(float)()
		self.mu_staco_id_theta = ROOT.std.vector(float)()
		self.mu_staco_id_qoverp = ROOT.std.vector(float)()

		self.mu_staco_id_theta_exPV = ROOT.std.vector(float)()
		self.mu_staco_id_qoverp_exPV = ROOT.std.vector(float)()

		self.mu_staco_trackd0pvunbiased = ROOT.std.vector(float)()
		self.mu_staco_trackz0pvunbiased = ROOT.std.vector(float)()
		self.mu_staco_tracksigd0pvunbiased = ROOT.std.vector(float)()
		self.mu_staco_tracksigz0pvunbiased = ROOT.std.vector(float)()

		self.mu_staco_ptcone20 = ROOT.std.vector(float)()
		self.mu_staco_ptcone30 = ROOT.std.vector(float)()
		self.mu_staco_ptcone40 = ROOT.std.vector(float)()

		self.mu_staco_etcone20 = ROOT.std.vector(float)()
		self.mu_staco_etcone30 = ROOT.std.vector(float)()
		self.mu_staco_etcone40 = ROOT.std.vector(float)()

		self.mu_staco_truth_type = ROOT.std.vector(int)()
		self.mu_staco_truth_mothertype = ROOT.std.vector(int)()
		self.mu_staco_truth_barcode = ROOT.std.vector(int)()
		self.mu_staco_truth_motherbarcode = ROOT.std.vector(int)()

		self.mu_staco_EFCB_dr = ROOT.std.vector(float)()
		self.mu_staco_EFCB_index = ROOT.std.vector(int)()

		##

		self.Tree1.Branch('mu_staco_n', self.mu_staco_n, 'mu_staco_n/I')

		self.Tree1.Branch('mu_staco_m', self.mu_staco_m)
		self.Tree1.Branch('mu_staco_E', self.mu_staco_E)
		self.Tree1.Branch('mu_staco_pt', self.mu_staco_pt)
		self.Tree1.Branch('mu_staco_eta', self.mu_staco_eta)
		self.Tree1.Branch('mu_staco_phi', self.mu_staco_phi)
		self.Tree1.Branch('mu_staco_charge', self.mu_staco_charge)
		self.Tree1.Branch('mu_staco_author', self.mu_staco_author)

		self.Tree1.Branch('mu_staco_loose', self.mu_staco_loose)
		self.Tree1.Branch('mu_staco_medium', self.mu_staco_medium)
		self.Tree1.Branch('mu_staco_tight', self.mu_staco_tight)

		self.Tree1.Branch('mu_staco_nBLHits', self.mu_staco_nBLHits)
		self.Tree1.Branch('mu_staco_nPixHits', self.mu_staco_nPixHits)
		self.Tree1.Branch('mu_staco_nSCTHits', self.mu_staco_nSCTHits)
		self.Tree1.Branch('mu_staco_nTRTHits', self.mu_staco_nTRTHits)
		self.Tree1.Branch('mu_staco_nTRTHighTHits', self.mu_staco_nTRTHighTHits)

		self.Tree1.Branch('mu_staco_nBLOutliers', self.mu_staco_nBLOutliers)
		self.Tree1.Branch('mu_staco_nPixOutliers', self.mu_staco_nPixOutliers)
		self.Tree1.Branch('mu_staco_nSCTOutliers', self.mu_staco_nSCTOutliers)
		self.Tree1.Branch('mu_staco_nTRTOutliers', self.mu_staco_nTRTOutliers)
		self.Tree1.Branch('mu_staco_nTRTHighTOutliers', self.mu_staco_nTRTHighTOutliers)

		self.Tree1.Branch('mu_staco_nPixHoles', self.mu_staco_nPixHoles)
		self.Tree1.Branch('mu_staco_nSCTHoles', self.mu_staco_nSCTHoles)
		self.Tree1.Branch('mu_staco_nTRTHoles', self.mu_staco_nTRTHoles)

		self.Tree1.Branch('mu_staco_expectBLayerHit', self.mu_staco_expectBLayerHit)

		self.Tree1.Branch('mu_staco_nPixDeadSensors', self.mu_staco_nPixDeadSensors)
		self.Tree1.Branch('mu_staco_nSCTDeadSensors', self.mu_staco_nSCTDeadSensors)

		self.Tree1.Branch('mu_staco_id_d0', self.mu_staco_id_d0)
		self.Tree1.Branch('mu_staco_id_z0', self.mu_staco_id_z0)
		self.Tree1.Branch('mu_staco_id_phi', self.mu_staco_id_phi)
		self.Tree1.Branch('mu_staco_id_theta', self.mu_staco_id_theta)
		self.Tree1.Branch('mu_staco_id_qoverp', self.mu_staco_id_qoverp)

		self.Tree1.Branch('mu_staco_id_theta_exPV', self.mu_staco_id_theta_exPV)
		self.Tree1.Branch('mu_staco_id_qoverp_exPV', self.mu_staco_id_qoverp_exPV)

		self.Tree1.Branch('mu_staco_trackd0pvunbiased', self.mu_staco_trackd0pvunbiased)
		self.Tree1.Branch('mu_staco_trackz0pvunbiased', self.mu_staco_trackz0pvunbiased)
		self.Tree1.Branch('mu_staco_tracksigd0pvunbiased', self.mu_staco_tracksigd0pvunbiased)
		self.Tree1.Branch('mu_staco_tracksigz0pvunbiased', self.mu_staco_tracksigz0pvunbiased)

		self.Tree1.Branch('mu_staco_ptcone20', self.mu_staco_ptcone20)
		self.Tree1.Branch('mu_staco_ptcone30', self.mu_staco_ptcone30)
		self.Tree1.Branch('mu_staco_ptcone40', self.mu_staco_ptcone40)

		self.Tree1.Branch('mu_staco_etcone20', self.mu_staco_etcone20)
		self.Tree1.Branch('mu_staco_etcone30', self.mu_staco_etcone30)
		self.Tree1.Branch('mu_staco_etcone40', self.mu_staco_etcone40)

		self.Tree1.Branch('mu_staco_truth_type', self.mu_staco_truth_type)
		self.Tree1.Branch('mu_staco_truth_mothertype', self.mu_staco_truth_mothertype)
		self.Tree1.Branch('mu_staco_truth_barcode', self.mu_staco_truth_barcode)
		self.Tree1.Branch('mu_staco_truth_motherbarcode', self.mu_staco_truth_motherbarcode)

		self.Tree1.Branch('mu_staco_EFCB_dr', self.mu_staco_EFCB_dr)
		self.Tree1.Branch('mu_staco_EFCB_index', self.mu_staco_EFCB_index)

		#########################
		# TRIGGERS		#
		#########################

		self.EF_e60_loose = array.array('i', [False])
		self.EF_e20_medium = array.array('i', [False])
		self.EF_e10_medium_mu6 = array.array('i', [False])
		self.EF_e10_medium_mu10 = array.array('i', [False])

		self.EF_2e10_medium = array.array('i', [False])
		self.EF_2e12_medium = array.array('i', [False])

		self.EF_2g15_loose = array.array('i', [False])
		self.EF_2g20_loose = array.array('i', [False])

		self.EF_mu18 = array.array('i', [False])
		self.EF_mu18_MG = array.array('i', [False])
		self.EF_mu20 = array.array('i', [False])
		self.EF_mu20_MG = array.array('i', [False])
		self.EF_mu22 = array.array('i', [False])
		self.EF_mu22_MG = array.array('i', [False])

		self.EF_mu20i = array.array('i', [False])

		self.EF_mu40_MSonly = array.array('i', [False])

		self.EF_2mu10_loose = array.array('i', [False])
		self.EF_mu15_mu10_EFFS = array.array('i', [False])

		##

		self.Tree2.Branch('EF_e60_loose', self.EF_e60_loose, 'EF_e60_loose/O')
		self.Tree2.Branch('EF_e20_medium', self.EF_e20_medium, 'EF_e20_medium/O')
		self.Tree2.Branch('EF_e10_medium_mu6', self.EF_e10_medium_mu6, 'EF_e10_medium_mu6/O')
		self.Tree2.Branch('EF_e10_medium_mu10', self.EF_e10_medium_mu10, 'EF_e10_medium_mu10/O')

		self.Tree2.Branch('EF_2e10_medium', self.EF_2e10_medium, 'EF_2e10_medium/O')
		self.Tree2.Branch('EF_2e12_medium', self.EF_2e12_medium, 'EF_2e12_medium/O')

		self.Tree2.Branch('EF_2g15_loose', self.EF_2g15_loose, 'EF_2g15_loose/O')
		self.Tree2.Branch('EF_2g20_loose', self.EF_2g20_loose, 'EF_2g20_loose/O')

		self.Tree2.Branch('EF_mu18', self.EF_mu18, 'EF_mu18/O')
		self.Tree2.Branch('EF_mu18_MG', self.EF_mu18_MG, 'EF_mu18_MG/O')
		self.Tree2.Branch('EF_mu20', self.EF_mu20, 'EF_mu20/O')
		self.Tree2.Branch('EF_mu20_MG', self.EF_mu20_MG, 'EF_mu20_MG/O')
		self.Tree2.Branch('EF_mu22', self.EF_mu22, 'EF_mu22/O')
		self.Tree2.Branch('EF_mu22_MG', self.EF_mu22_MG, 'EF_mu22_MG/O')

		self.Tree2.Branch('EF_mu20i', self.EF_mu20i, 'EF_mu20i/O')

		self.Tree2.Branch('EF_mu40_MSonly', self.EF_mu40_MSonly, 'EF_mu40_MSonly/O')

		self.Tree2.Branch('EF_2mu10_loose', self.EF_2mu10_loose, 'EF_2mu10_loose/O')
		self.Tree2.Branch('EF_mu15_mu10_EFFS', self.EF_mu15_mu10_EFFS, 'EF_mu15_mu10_EFFS/O')

		#########################
		# TRIGGER ELECTRONS	#
		#########################

		self.trig_EF_el_n = array.array('i', [0])

		self.trig_EF_el_eta = ROOT.std.vector(float)()
		self.trig_EF_el_phi = ROOT.std.vector(float)()

		self.trig_EF_el_EF_e60_loose = ROOT.std.vector(int)()
		self.trig_EF_el_EF_e20_medium = ROOT.std.vector(int)()
		self.trig_EF_el_EF_e10_medium_mu6 = ROOT.std.vector(int)()
		self.trig_EF_el_EF_e10_medium_mu10 = ROOT.std.vector(int)()

		self.trig_EF_el_EF_2e10_medium = ROOT.std.vector(int)()
		self.trig_EF_el_EF_2e12_medium = ROOT.std.vector(int)()

		self.trig_EF_el_EF_2g15_loose = ROOT.std.vector(int)()
		self.trig_EF_el_EF_2g20_loose = ROOT.std.vector(int)()

		##

		self.Tree2.Branch('trig_EF_el_n', self.trig_EF_el_n, 'trig_EF_el_n/I')

		self.Tree2.Branch('trig_EF_el_eta', self.trig_EF_el_eta)
		self.Tree2.Branch('trig_EF_el_phi', self.trig_EF_el_phi)

		self.Tree2.Branch('trig_EF_el_EF_e60_loose', self.trig_EF_el_EF_e60_loose)
		self.Tree2.Branch('trig_EF_el_EF_e20_medium', self.trig_EF_el_EF_e20_medium)
		self.Tree2.Branch('trig_EF_el_EF_e10_medium_mu6', self.trig_EF_el_EF_e10_medium_mu6)
		self.Tree2.Branch('trig_EF_el_EF_e10_medium_mu10', self.trig_EF_el_EF_e10_medium_mu10)

		self.Tree2.Branch('trig_EF_el_EF_2e10_medium', self.trig_EF_el_EF_2e10_medium)
		self.Tree2.Branch('trig_EF_el_EF_2e12_medium', self.trig_EF_el_EF_2e12_medium)

		self.Tree2.Branch('trig_EF_el_EF_2g15_loose', self.trig_EF_el_EF_2g15_loose)
		self.Tree2.Branch('trig_EF_el_EF_2g20_loose', self.trig_EF_el_EF_2g20_loose)

		#########################
		# TRIGGER MUONS		#
		#########################

		self.trig_EF_trigmuonef_n = array.array('i', [0])

		self.trig_EF_trigmuonef_EF_mu18 = ROOT.std.vector(int)()
		self.trig_EF_trigmuonef_EF_mu18_MG = ROOT.std.vector(int)()
		self.trig_EF_trigmuonef_EF_mu20 = ROOT.std.vector(int)()
		self.trig_EF_trigmuonef_EF_mu20_MG = ROOT.std.vector(int)()
		self.trig_EF_trigmuonef_EF_mu22 = ROOT.std.vector(int)()
		self.trig_EF_trigmuonef_EF_mu22_MG = ROOT.std.vector(int)()

		self.trig_EF_trigmuonef_EF_mu20i = ROOT.std.vector(int)()

		self.trig_EF_trigmuonef_EF_2mu10_loose = ROOT.std.vector(int)()
		self.trig_EF_trigmuonef_EF_mu15_mu10_EFFS = ROOT.std.vector(int)()

		##

		self.Tree2.Branch('trig_EF_trigmuonef_n', self.trig_EF_trigmuonef_n, 'trig_EF_trigmuonef_n/I')

		self.Tree2.Branch('trig_EF_trigmuonef_EF_mu18', self.trig_EF_trigmuonef_EF_mu18)
		self.Tree2.Branch('trig_EF_trigmuonef_EF_mu18_MG', self.trig_EF_trigmuonef_EF_mu18_MG)
		self.Tree2.Branch('trig_EF_trigmuonef_EF_mu20', self.trig_EF_trigmuonef_EF_mu20)
		self.Tree2.Branch('trig_EF_trigmuonef_EF_mu20_MG', self.trig_EF_trigmuonef_EF_mu20_MG)
		self.Tree2.Branch('trig_EF_trigmuonef_EF_mu22', self.trig_EF_trigmuonef_EF_mu22)
		self.Tree2.Branch('trig_EF_trigmuonef_EF_mu22_MG', self.trig_EF_trigmuonef_EF_mu22_MG)

		self.Tree2.Branch('trig_EF_trigmuonef_EF_mu20i', self.trig_EF_trigmuonef_EF_mu20i)

		self.Tree2.Branch('trig_EF_trigmuonef_EF_2mu10_loose', self.trig_EF_trigmuonef_EF_2mu10_loose)
		self.Tree2.Branch('trig_EF_trigmuonef_EF_mu15_mu10_EFFS', self.trig_EF_trigmuonef_EF_mu15_mu10_EFFS)

	#####################################################################

	def treeCleaner(self):
		#########################
		# EVENT			#
		#########################

		self.RunNumber[0] = 0
		self.EventNumber[0] = 0
		self.lbn[0] = 0

		self.pixelError[0] = 0
		self.sctError[0] = 0
		self.trtError[0] = 0
		self.larError[0] = 0
		self.muonError[0] = 0

		self.mcevt_weight.clear()

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

		self.el_m.clear()
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

		self.el_E237.clear()
		self.el_E277.clear()
		self.el_emaxs1.clear()
		self.el_Emax2.clear()
		self.el_Ethad.clear()
		self.el_Ethad1.clear()
		self.el_f1.clear()
		self.el_weta2.clear()
		self.el_wstot.clear()
		self.el_etap.clear()

		self.el_nBLHits.clear()
		self.el_nPixHits.clear()
		self.el_nSCTHits.clear()
		self.el_nTRTHits.clear()
		self.el_nTRTHighTHits.clear()

		self.el_nBLOutliers.clear()
		self.el_nPixOutliers.clear()
		self.el_nSCTOutliers.clear()
		self.el_nTRTOutliers.clear()
		self.el_nTRTHighTOutliers.clear()

		self.el_nPixHoles.clear()
		self.el_nSCTHoles.clear()
		self.el_nTRTHoles.clear()

		self.el_expectBLayerHit.clear()

		self.el_deltaeta1.clear()
		self.el_deltaeta2.clear()

		self.el_trackd0.clear()
		self.el_trackz0.clear()
		self.el_trackpt.clear()
		self.el_tracketa.clear()
		self.el_trackphi.clear()
		self.el_tracktheta.clear()
		self.el_trackqoverp.clear()

		self.el_trackd0pvunbiased.clear()
		self.el_trackz0pvunbiased.clear()
		self.el_tracksigd0pvunbiased.clear()
		self.el_tracksigz0pvunbiased.clear()

		self.el_ptcone20.clear()
		self.el_ptcone30.clear()
		self.el_ptcone40.clear()

		self.el_Etcone20.clear()
		self.el_Etcone30.clear()
		self.el_Etcone40.clear()

		self.el_truth_type.clear()
		self.el_truth_barcode.clear()
		self.el_truth_mothertype.clear()
		self.el_truth_motherbarcode.clear()

		self.el_type.clear()
		self.el_origin.clear()
		self.el_typebkg.clear()
		self.el_originbkg.clear()

		self.el_EF_dr.clear()
		self.el_EF_index.clear()

		#########################
		# MUONS MUID		#
		#########################

		self.mu_muid_n[0] = 0

		self.mu_muid_m.clear()
		self.mu_muid_E.clear()
		self.mu_muid_pt.clear()
		self.mu_muid_eta.clear()
		self.mu_muid_phi.clear()
		self.mu_muid_charge.clear()
		self.mu_muid_author.clear()

		self.mu_muid_loose.clear()
		self.mu_muid_medium.clear()
		self.mu_muid_tight.clear()

		self.mu_muid_nBLHits.clear()
		self.mu_muid_nPixHits.clear()
		self.mu_muid_nSCTHits.clear()
		self.mu_muid_nTRTHits.clear()
		self.mu_muid_nTRTHighTHits.clear()

		self.mu_muid_nBLOutliers.clear()
		self.mu_muid_nPixOutliers.clear()
		self.mu_muid_nSCTOutliers.clear()
		self.mu_muid_nTRTOutliers.clear()
		self.mu_muid_nTRTHighTOutliers.clear()

		self.mu_muid_nPixHoles.clear()
		self.mu_muid_nSCTHoles.clear()
		self.mu_muid_nTRTHoles.clear()

		self.mu_muid_expectBLayerHit.clear()

		self.mu_muid_nPixDeadSensors.clear()
		self.mu_muid_nSCTDeadSensors.clear()

		self.mu_muid_id_d0.clear()
		self.mu_muid_id_z0.clear()
		self.mu_muid_id_phi.clear()
		self.mu_muid_id_theta.clear()
		self.mu_muid_id_qoverp.clear()

		self.mu_muid_id_theta_exPV.clear()
		self.mu_muid_id_qoverp_exPV.clear()

		self.mu_muid_trackd0pvunbiased.clear()
		self.mu_muid_trackz0pvunbiased.clear()
		self.mu_muid_tracksigd0pvunbiased.clear()
		self.mu_muid_tracksigz0pvunbiased.clear()

		self.mu_muid_ptcone20.clear()
		self.mu_muid_ptcone30.clear()
		self.mu_muid_ptcone40.clear()

		self.mu_muid_etcone20.clear()
		self.mu_muid_etcone30.clear()
		self.mu_muid_etcone40.clear()

		self.mu_muid_truth_type.clear()
		self.mu_muid_truth_barcode.clear()
		self.mu_muid_truth_mothertype.clear()
		self.mu_muid_truth_motherbarcode.clear()

		self.mu_muid_EFCB_dr.clear()
		self.mu_muid_EFCB_index.clear()

		#########################
		# MUONS STACO		#
		#########################

		self.mu_staco_n[0] = 0

		self.mu_staco_m.clear()
		self.mu_staco_E.clear()
		self.mu_staco_pt.clear()
		self.mu_staco_eta.clear()
		self.mu_staco_phi.clear()
		self.mu_staco_charge.clear()
		self.mu_staco_author.clear()

		self.mu_staco_loose.clear()
		self.mu_staco_medium.clear()
		self.mu_staco_tight.clear()

		self.mu_staco_nBLHits.clear()
		self.mu_staco_nPixHits.clear()
		self.mu_staco_nSCTHits.clear()
		self.mu_staco_nTRTHits.clear()
		self.mu_staco_nTRTHighTHits.clear()

		self.mu_staco_nBLOutliers.clear()
		self.mu_staco_nPixOutliers.clear()
		self.mu_staco_nSCTOutliers.clear()
		self.mu_staco_nTRTOutliers.clear()
		self.mu_staco_nTRTHighTOutliers.clear()

		self.mu_staco_nPixHoles.clear()
		self.mu_staco_nSCTHoles.clear()
		self.mu_staco_nTRTHoles.clear()

		self.mu_staco_expectBLayerHit.clear()

		self.mu_staco_nPixDeadSensors.clear()
		self.mu_staco_nSCTDeadSensors.clear()

		self.mu_staco_id_d0.clear()
		self.mu_staco_id_z0.clear()
		self.mu_staco_id_phi.clear()
		self.mu_staco_id_theta.clear()
		self.mu_staco_id_qoverp.clear()

		self.mu_staco_id_theta_exPV.clear()
		self.mu_staco_id_qoverp_exPV.clear()

		self.mu_staco_trackd0pvunbiased.clear()
		self.mu_staco_trackz0pvunbiased.clear()
		self.mu_staco_tracksigd0pvunbiased.clear()
		self.mu_staco_tracksigz0pvunbiased.clear()

		self.mu_staco_ptcone20.clear()
		self.mu_staco_ptcone30.clear()
		self.mu_staco_ptcone40.clear()

		self.mu_staco_etcone20.clear()
		self.mu_staco_etcone30.clear()
		self.mu_staco_etcone40.clear()

		self.mu_staco_truth_type.clear()
		self.mu_staco_truth_barcode.clear()
		self.mu_staco_truth_mothertype.clear()
		self.mu_staco_truth_motherbarcode.clear()

		self.mu_staco_EFCB_dr.clear()
		self.mu_staco_EFCB_index.clear()

		#########################
		# TRIGGERS		#
		#########################

		self.EF_e60_loose[0] = False
		self.EF_e20_medium[0] = False
		self.EF_e10_medium_mu6[0] = False
		self.EF_e10_medium_mu10[0] = False

		self.EF_2e10_medium[0] = False
		self.EF_2e12_medium[0] = False

		self.EF_2g15_loose[0] = False
		self.EF_2g20_loose[0] = False

		##

		self.EF_mu18[0] = False
		self.EF_mu18_MG[0] = False
		self.EF_mu20[0] = False
		self.EF_mu20_MG[0] = False
		self.EF_mu22[0] = False
		self.EF_mu22_MG[0] = False

		self.EF_mu20i[0] = False

		self.EF_mu40_MSonly[0] = False

		self.EF_2mu10_loose[0] = False
		self.EF_mu15_mu10_EFFS[0] = False

		#########################
		# TRIGGER ELECTRONS	#
		#########################

		self.trig_EF_el_n[0] = 0

		self.trig_EF_el_eta.clear()
		self.trig_EF_el_phi.clear()

		self.trig_EF_el_EF_e60_loose.clear()
		self.trig_EF_el_EF_e20_medium.clear()
		self.trig_EF_el_EF_e10_medium_mu6.clear()
		self.trig_EF_el_EF_e10_medium_mu10.clear()

		self.trig_EF_el_EF_2e10_medium.clear()
		self.trig_EF_el_EF_2e12_medium.clear()

		self.trig_EF_el_EF_2g15_loose.clear()
		self.trig_EF_el_EF_2g20_loose.clear()

		#########################
		# TRIGGER MUONS		#
		#########################

		self.trig_EF_trigmuonef_n[0] = 0

		self.trig_EF_trigmuonef_EF_mu18.clear()
		self.trig_EF_trigmuonef_EF_mu18_MG.clear()
		self.trig_EF_trigmuonef_EF_mu20.clear()
		self.trig_EF_trigmuonef_EF_mu20_MG.clear()
		self.trig_EF_trigmuonef_EF_mu22.clear()
		self.trig_EF_trigmuonef_EF_mu22_MG.clear()

		self.trig_EF_trigmuonef_EF_mu20i.clear()

		self.trig_EF_trigmuonef_EF_2mu10_loose.clear()
		self.trig_EF_trigmuonef_EF_mu15_mu10_EFFS.clear()

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
		elif 'ByteStreamEventInfo' in L:
			event = self.StoreGateSvc['ByteStreamEventInfo']
		else:
			return PyAthena.StatusCode.Failure

		eventID = event.event_ID()
		eventType = event.event_type()

		#############################################################

		self.RunNumber[0] = eventID.run_number()
		self.EventNumber[0] = eventID.event_number()
		self.lbn[0] = eventID.lumi_block()

		self.pixelError[0] = event.errorState(PyAthena.EventInfo.Pixel)
		self.sctError[0] = event.errorState(PyAthena.EventInfo.SCT)
		self.trtError[0] = event.errorState(PyAthena.EventInfo.TRT)
		self.larError[0] = event.errorState(PyAthena.EventInfo.LAr)
		self.muonError[0] = event.errorState(PyAthena.EventInfo.Muon)

		self.mcevt_weight.push_back(eventType.mc_event_weight())

		#############################################################
		# PRIMARY VERTICES & PRIMARY TRACKS			    #
		#############################################################

		vertices = self.StoreGateSvc['VxPrimaryCandidate']

		#############################################################

		NPV = 0

		for vertex in vertices:
			self.vxp_x.push_back(vertex.recVertex().position().x())
			self.vxp_y.push_back(vertex.recVertex().position().y())
			self.vxp_z.push_back(vertex.recVertex().position().z())

			self.vxp_nTracks.push_back(len(vertex.vxTrackAtVertex()))

			##

			NPV += 1
			self.vxp_n[0] += 1

		#############################################################
		# SKIMMING						    #
		#############################################################

		N1 = 0
		N2 = 0
		N3 = 0

		for obj in self.StoreGateSvc['ElectronAODCollection']:
			if (obj.author() == 1 or obj.author() == 3) and obj.isElectron(PyAthena.egammaPID.ElectronMedium):
				N1 = N1 + 1

		for obj in self.StoreGateSvc['MuidMuonCollection']:
			if obj.isTight():
				N2 = N2 + 1

		for obj in self.StoreGateSvc['StacoMuonCollection']:
			if obj.author() == 6 or obj.author() == 7:
				N3 = N3 + 1

		#############################################################

		if NPV == 0 or (N1 < 2 and N2 < 2 and N3 < 2):
			return PyAthena.StatusCode.Success

		#############################################################
		# ELECTRONS						    #
		#############################################################

		electrons = self.StoreGateSvc['ElectronAODCollection']

		#############################################################

		for electron in electrons:

			m = electron.m()
			E = electron.e()
			Et = electron.et()
			pt = electron.pt()
			eta = electron.eta()
			phi = electron.phi()
			charge = electron.charge()
			author = electron.author()
			isEM = electron.isem()

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
				trackpt = track.pt()
				tracketa = track.eta()

				summary = track.trackSummary()

				if summary:
					nBLHits = summary.get(ROOT.Trk.numberOfBLayerHits)
					nPixHits = summary.get(ROOT.Trk.numberOfPixelHits)
					nSCTHits = summary.get(ROOT.Trk.numberOfSCTHits)
					nTRTHits = summary.get(ROOT.Trk.numberOfTRTHits)
					nTRTHighTHits = summary.get(ROOT.Trk.numberOfTRTHighThresholdHits)

					nBLOutliers = summary.get(ROOT.Trk.numberOfBLayerOutliers)
					nPixOutliers = summary.get(ROOT.Trk.numberOfPixelOutliers)
					nSCTOutliers = summary.get(ROOT.Trk.numberOfSCTOutliers)
					nTRTOutliers = summary.get(ROOT.Trk.numberOfTRTOutliers)
					nTRTHighTOutliers = summary.get(ROOT.Trk.numberOfTRTHighThresholdOutliers)

					nPixHoles = summary.get(ROOT.Trk.numberOfPixelHoles)
					nSCTHoles = summary.get(ROOT.Trk.numberOfSCTHoles)
					nTRTHoles = summary.get(ROOT.Trk.numberOfTRTHoles)
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

			if electron.nDetails() != 0:
				E237 = electron.detailValue(ROOT.egammaParameters.e237)
				E277 = electron.detailValue(ROOT.egammaParameters.e277)
				emaxs1 = electron.detailValue(ROOT.egammaParameters.emaxs1)
				Emax2 = electron.detailValue(ROOT.egammaParameters.e2tsts1)
				Ethad = electron.detailValue(ROOT.egammaParameters.ethad)
				Ethad1 = electron.detailValue(ROOT.egammaParameters.ethad1)
				f1 = electron.detailValue(ROOT.egammaParameters.f1)
				weta2 = electron.detailValue(ROOT.egammaParameters.weta2)
				wstot = electron.detailValue(ROOT.egammaParameters.wtots1)
				etap = electron.detailValue(ROOT.egammaParameters.etap)

				ptcone20 = electron.detailValue(ROOT.egammaParameters.ptcone20)
				ptcone30 = electron.detailValue(ROOT.egammaParameters.ptcone30)
				ptcone40 = electron.detailValue(ROOT.egammaParameters.ptcone40)

				Etcone20 = electron.detailValue(ROOT.egammaParameters.etcone20)
				Etcone30 = electron.detailValue(ROOT.egammaParameters.etcone30)
				Etcone40 = electron.detailValue(ROOT.egammaParameters.etcone40)

				expectBLayerHit = electron.detailValue(ROOT.egammaParameters.expectHitInBLayer) > 0

				deltaeta1 = electron.detailValue(ROOT.egammaParameters.deltaEta1)
				deltaeta2 = electron.detailValue(ROOT.egammaParameters.deltaEta2)
			else:
				continue

			##

			if isMC == False:
				truth_type = 0
				truth_mothertype = 0
				truth_barcode = 0
				truth_motherbarcode = 0

				epyt = 0
				origin = 0
				typebkg = 0
				originbkg = 0
			else:
				truth = self.StoreGateSvc['SpclMC']

				index = particleMatching(eta, phi, truth)[0]

				if index >= 0:
					truth_type = truth[index].pdgId()
					truth_barcode = truth[index].barcode()

					if truth[index].mother():
						truth_mothertype = truth[index].mother().pdgId()
						truth_motherbarcode = truth[index].mother().barcode()
					else:
						truth_mothertype = 0
						truth_motherbarcode = 0
				else:
					truth_type = 0
					truth_barcode = 0
					truth_mothertype = 0
					truth_motherbarcode = 0

				##

				epyt, origin = self.MCTruthClassifier.particleTruthClassifier(electron)

				if epyt > 0\
				   and     \
				   epyt < 5:
					mc_electron = self.MCTruthClassifier.getGenPart()

					typebkg, originbkg = self.MCTruthClassifier.checkOrigOfBkgElec(mc_electron)
				else:
					typebkg, originbkg = 0x000000000000000000000000, 0x000000000000000000000000

				#####################################
				# PATCH FOR JF17		    #
				#####################################

				if self.RunNumber == 105802 and ((epyt == 2 and origin == 13) or (epyt == 4 and (originbkg == 13 or originbkg == 40))): continue

				#####################################
				# PATCH FOR JF17		    #
				#####################################

			##

			if (nPixHits + nSCTHits) >= 4:
				p = particleMatching(tracketa, trackphi, self.StoreGateSvc['HLT_egamma_Electrons'])
			else:
				p = particleMatching( cl_eta ,  cl_phi , self.StoreGateSvc['HLT_egamma_Electrons'])

			EF_index = p[0]
			EF_dr = p[1]

			##

			self.el_m.push_back(m)
			self.el_E.push_back(E)
			self.el_Et.push_back(Et)
			self.el_pt.push_back(pt)
			self.el_eta.push_back(eta)
			self.el_phi.push_back(phi)
			self.el_charge.push_back(charge)
			self.el_author.push_back(author)
			self.el_isEM.push_back(isEM)
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

			self.el_E237.push_back(E237)
			self.el_E277.push_back(E277)
			self.el_emaxs1.push_back(emaxs1)
			self.el_Emax2.push_back(Emax2)
			self.el_Ethad.push_back(Ethad)
			self.el_Ethad1.push_back(Ethad1)
			self.el_f1.push_back(f1)
			self.el_weta2.push_back(weta2)
			self.el_wstot.push_back(wstot)
			self.el_etap.push_back(etap)

			self.el_nBLHits.push_back(nBLHits)
			self.el_nPixHits.push_back(nPixHits)
			self.el_nSCTHits.push_back(nSCTHits)
			self.el_nTRTHits.push_back(nTRTHits)
			self.el_nTRTHighTHits.push_back(nTRTHighTHits)

			self.el_nBLOutliers.push_back(nBLOutliers)
			self.el_nPixOutliers.push_back(nPixOutliers)
			self.el_nSCTOutliers.push_back(nSCTOutliers)
			self.el_nTRTOutliers.push_back(nTRTOutliers)
			self.el_nTRTHighTOutliers.push_back(nTRTHighTOutliers)

			self.el_nPixHoles.push_back(nPixHoles)
			self.el_nSCTHoles.push_back(nSCTHoles)
			self.el_nTRTHoles.push_back(nTRTHoles)

			self.el_expectBLayerHit.push_back(expectBLayerHit)

			self.el_deltaeta1.push_back(deltaeta1)
			self.el_deltaeta2.push_back(deltaeta2)

			self.el_trackd0.push_back(trackd0)
			self.el_trackz0.push_back(trackz0)
			self.el_trackpt.push_back(trackpt)
			self.el_tracketa.push_back(tracketa)
			self.el_trackphi.push_back(trackphi)
			self.el_tracktheta.push_back(tracktheta)
			self.el_trackqoverp.push_back(trackqoverp)

			self.el_trackd0pvunbiased.push_back(trackd0pvunbiased)
			self.el_trackz0pvunbiased.push_back(trackz0pvunbiased)
			self.el_tracksigd0pvunbiased.push_back(tracksigd0pvunbiased)
			self.el_tracksigz0pvunbiased.push_back(tracksigz0pvunbiased)

			self.el_ptcone20.push_back(ptcone20)
			self.el_ptcone30.push_back(ptcone30)
			self.el_ptcone40.push_back(ptcone40)

			self.el_Etcone20.push_back(Etcone20)
			self.el_Etcone30.push_back(Etcone30)
			self.el_Etcone40.push_back(Etcone40)

			self.el_truth_type.push_back(truth_type)
			self.el_truth_barcode.push_back(truth_barcode)
			self.el_truth_mothertype.push_back(truth_mothertype)
			self.el_truth_motherbarcode.push_back(truth_motherbarcode)

			self.el_type.push_back(epyt)
			self.el_origin.push_back(origin)
			self.el_typebkg.push_back(typebkg)
			self.el_originbkg.push_back(originbkg)

			self.el_EF_dr.push_back(EF_dr)
			self.el_EF_index.push_back(EF_index)

			##

			self.el_n[0] += 1

		#############################################################
		# MUONS MUID						    #
		#############################################################

		muons = self.StoreGateSvc['MuidMuonCollection']

		#############################################################

		for muon in muons:

			m = muon.m()
			E = muon.e()
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

			track = muon.track()

			if track:
				summary = track.trackSummary()

				if summary:
					nBLHits = summary.get(ROOT.Trk.numberOfBLayerHits)
					nPixHits = summary.get(ROOT.Trk.numberOfPixelHits)
					nSCTHits = summary.get(ROOT.Trk.numberOfSCTHits)
					nTRTHits = summary.get(ROOT.Trk.numberOfTRTHits)
					nTRTHighTHits = summary.get(ROOT.Trk.numberOfTRTHighThresholdHits)

					nBLOutliers = summary.get(ROOT.Trk.numberOfBLayerOutliers)
					nPixOutliers = summary.get(ROOT.Trk.numberOfPixelOutliers)
					nSCTOutliers = summary.get(ROOT.Trk.numberOfSCTOutliers)
					nTRTOutliers = summary.get(ROOT.Trk.numberOfTRTOutliers)
					nTRTHighTOutliers = summary.get(ROOT.Trk.numberOfTRTHighThresholdOutliers)

					nPixHoles = summary.get(ROOT.Trk.numberOfPixelHoles)
					nSCTHoles = summary.get(ROOT.Trk.numberOfSCTHoles)
					nTRTHoles = summary.get(ROOT.Trk.numberOfTRTHoles)

					expectBLayerHit = summary.get(ROOT.Trk.expectBLayerHit)

					nPixDeadSensors = summary.get(ROOT.Trk.numberOfPixelDeadSensors)
					nSCTDeadSensors = summary.get(ROOT.Trk.numberOfSCTDeadSensors)
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

			track = muon.inDetTrackParticle()

			if track:
				perigee = track.measuredPerigee()

				if perigee:
					id_d0 = perigee.parameters()[0]
					id_z0 = perigee.parameters()[1]
					id_phi = perigee.parameters()[2]
					id_theta = perigee.parameters()[3]
					id_qoverp = perigee.parameters()[4]
				else:
					continue
			else:
				continue

			##

			ptcone20 = muon.parameter(PyAthena.MuonParameters.ptcone20)
			ptcone30 = muon.parameter(PyAthena.MuonParameters.ptcone30)
			ptcone40 = muon.parameter(PyAthena.MuonParameters.ptcone40)

			Etcone20 = muon.parameter(PyAthena.MuonParameters.etcone20)
			Etcone30 = muon.parameter(PyAthena.MuonParameters.etcone30)
			Etcone40 = muon.parameter(PyAthena.MuonParameters.etcone40)

			##

			if isMC == False:
				truth_type = 0
				truth_mothertype = 0
				truth_barcode = 0
				truth_motherbarcode = 0
			else:
				truth = self.StoreGateSvc['SpclMC']

				index = particleMatching(eta, phi, truth)[0]

				if index >= 0:
					truth_type = truth[index].pdgId()
					truth_barcode = truth[index].barcode()

					if truth[index].mother():
						truth_mothertype = truth[index].mother().pdgId()
						truth_motherbarcode = truth[index].mother().barcode()
					else:
						truth_mothertype = 0
						truth_motherbarcode = 0
				else:
					truth_type = 0
					truth_barcode = 0
					truth_mothertype = 0
					truth_motherbarcode = 0

			##

			p = particleMatching(eta, phi, self.StoreGateSvc['HLT_MuonEFInfo'])
			EF_index = p[0]
			EF_dr = p[1]

			##

			self.mu_muid_m.push_back(m)
			self.mu_muid_E.push_back(E)
			self.mu_muid_pt.push_back(pt)
			self.mu_muid_eta.push_back(eta)
			self.mu_muid_phi.push_back(phi)
			self.mu_muid_charge.push_back(charge)
			self.mu_muid_author.push_back(author)

			self.mu_muid_loose.push_back(loose)
			self.mu_muid_medium.push_back(medium)
			self.mu_muid_tight.push_back(tight)

			self.mu_muid_nBLHits.push_back(nBLHits)
			self.mu_muid_nPixHits.push_back(nPixHits)
			self.mu_muid_nSCTHits.push_back(nSCTHits)
			self.mu_muid_nTRTHits.push_back(nTRTHits)
			self.mu_muid_nTRTHighTHits.push_back(nTRTHighTHits)

			self.mu_muid_nBLOutliers.push_back(nBLOutliers)
			self.mu_muid_nPixOutliers.push_back(nPixOutliers)
			self.mu_muid_nSCTOutliers.push_back(nSCTOutliers)
			self.mu_muid_nTRTOutliers.push_back(nTRTOutliers)
			self.mu_muid_nTRTHighTOutliers.push_back(nTRTHighTOutliers)

			self.mu_muid_nPixHoles.push_back(nPixHoles)
			self.mu_muid_nSCTHoles.push_back(nSCTHoles)
			self.mu_muid_nTRTHoles.push_back(nTRTHoles)

			self.mu_muid_expectBLayerHit.push_back(expectBLayerHit)

			self.mu_muid_nPixDeadSensors.push_back(nPixDeadSensors)
			self.mu_muid_nSCTDeadSensors.push_back(nSCTDeadSensors)

			self.mu_muid_id_d0.push_back(id_d0)
			self.mu_muid_id_z0.push_back(id_z0)
			self.mu_muid_id_phi.push_back(id_phi)
			self.mu_muid_id_theta.push_back(id_theta)
			self.mu_muid_id_qoverp.push_back(id_qoverp)

			self.mu_muid_id_theta_exPV.push_back(id_theta)
			self.mu_muid_id_qoverp_exPV.push_back(id_qoverp)

			self.mu_muid_trackd0pvunbiased.push_back(trackd0pvunbiased)
			self.mu_muid_trackz0pvunbiased.push_back(trackz0pvunbiased)
			self.mu_muid_tracksigd0pvunbiased.push_back(tracksigd0pvunbiased)
			self.mu_muid_tracksigz0pvunbiased.push_back(tracksigz0pvunbiased)

			self.mu_muid_ptcone20.push_back(ptcone20)
			self.mu_muid_ptcone30.push_back(ptcone30)
			self.mu_muid_ptcone40.push_back(ptcone40)

			self.mu_muid_etcone20.push_back(Etcone20)
			self.mu_muid_etcone30.push_back(Etcone30)
			self.mu_muid_etcone40.push_back(Etcone40)

			self.mu_muid_truth_type.push_back(truth_type)
			self.mu_muid_truth_barcode.push_back(truth_barcode)
			self.mu_muid_truth_mothertype.push_back(truth_mothertype)
			self.mu_muid_truth_motherbarcode.push_back(truth_motherbarcode)

			self.mu_muid_EFCB_dr.push_back(EF_dr)
			self.mu_muid_EFCB_index.push_back(EF_index)

			##

			self.mu_muid_n[0] += 1

		#############################################################
		# MUONS STACO						    #
		#############################################################

		muons = self.StoreGateSvc['StacoMuonCollection']

		#############################################################

		for muon in muons:

			m = muon.m()
			E = muon.e()
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

			track = muon.track()

			if track:
				summary = track.trackSummary()

				if summary:
					nBLHits = summary.get(ROOT.Trk.numberOfBLayerHits)
					nPixHits = summary.get(ROOT.Trk.numberOfPixelHits)
					nSCTHits = summary.get(ROOT.Trk.numberOfSCTHits)
					nTRTHits = summary.get(ROOT.Trk.numberOfTRTHits)
					nTRTHighTHits = summary.get(ROOT.Trk.numberOfTRTHighThresholdHits)

					nBLOutliers = summary.get(ROOT.Trk.numberOfBLayerOutliers)
					nPixOutliers = summary.get(ROOT.Trk.numberOfPixelOutliers)
					nSCTOutliers = summary.get(ROOT.Trk.numberOfSCTOutliers)
					nTRTOutliers = summary.get(ROOT.Trk.numberOfTRTOutliers)
					nTRTHighTOutliers = summary.get(ROOT.Trk.numberOfTRTHighThresholdOutliers)

					nPixHoles = summary.get(ROOT.Trk.numberOfPixelHoles)
					nSCTHoles = summary.get(ROOT.Trk.numberOfSCTHoles)
					nTRTHoles = summary.get(ROOT.Trk.numberOfTRTHoles)

					expectBLayerHit = summary.get(ROOT.Trk.expectBLayerHit)

					nPixDeadSensors = summary.get(ROOT.Trk.numberOfPixelDeadSensors)
					nSCTDeadSensors = summary.get(ROOT.Trk.numberOfSCTDeadSensors)
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

			track = muon.inDetTrackParticle()

			if track:
				perigee = track.measuredPerigee()

				if perigee:
					trackd0 = perigee.parameters()[0]
					trackz0 = perigee.parameters()[1]
					trackphi = perigee.parameters()[2]
					tracktheta = perigee.parameters()[3]
					trackqoverp = perigee.parameters()[4]
				else:
					continue
			else:
				continue

			##

			ptcone20 = muon.parameter(PyAthena.MuonParameters.ptcone20)
			ptcone30 = muon.parameter(PyAthena.MuonParameters.ptcone30)
			ptcone40 = muon.parameter(PyAthena.MuonParameters.ptcone40)

			Etcone20 = muon.parameter(PyAthena.MuonParameters.etcone20)
			Etcone30 = muon.parameter(PyAthena.MuonParameters.etcone30)
			Etcone40 = muon.parameter(PyAthena.MuonParameters.etcone40)

			##

			if isMC == False:
				truth_type = 0
				truth_mothertype = 0
				truth_barcode = 0
				truth_motherbarcode = 0
			else:
				truth = self.StoreGateSvc['SpclMC']

				index = particleMatching(eta, phi, truth)[0]

				if index >= 0:
					truth_type = truth[index].pdgId()
					truth_barcode = truth[index].barcode()

					if truth[index].mother():
						truth_mothertype = truth[index].mother().pdgId()
						truth_motherbarcode = truth[index].mother().barcode()
					else:
						truth_mothertype = 0
						truth_motherbarcode = 0
				else:
					truth_type = 0
					truth_barcode = 0
					truth_mothertype = 0
					truth_motherbarcode = 0

			##

			p = particleMatching(eta, phi, self.StoreGateSvc['HLT_MuonEFInfo'])
			EF_index = p[0]
			EF_dr = p[1]

			##

			self.mu_staco_m.push_back(m)
			self.mu_staco_E.push_back(E)
			self.mu_staco_pt.push_back(pt)
			self.mu_staco_eta.push_back(eta)
			self.mu_staco_phi.push_back(phi)
			self.mu_staco_charge.push_back(charge)
			self.mu_staco_author.push_back(author)

			self.mu_staco_loose.push_back(loose)
			self.mu_staco_medium.push_back(medium)
			self.mu_staco_tight.push_back(tight)

			self.mu_staco_nBLHits.push_back(nBLHits)
			self.mu_staco_nPixHits.push_back(nPixHits)
			self.mu_staco_nSCTHits.push_back(nSCTHits)
			self.mu_staco_nTRTHits.push_back(nTRTHits)
			self.mu_staco_nTRTHighTHits.push_back(nTRTHighTHits)

			self.mu_staco_nBLOutliers.push_back(nBLOutliers)
			self.mu_staco_nPixOutliers.push_back(nPixOutliers)
			self.mu_staco_nSCTOutliers.push_back(nSCTOutliers)
			self.mu_staco_nTRTOutliers.push_back(nTRTOutliers)
			self.mu_staco_nTRTHighTOutliers.push_back(nTRTHighTOutliers)

			self.mu_staco_nPixHoles.push_back(nPixHoles)
			self.mu_staco_nSCTHoles.push_back(nSCTHoles)
			self.mu_staco_nTRTHoles.push_back(nTRTHoles)

			self.mu_staco_expectBLayerHit.push_back(expectBLayerHit)

			self.mu_staco_nPixDeadSensors.push_back(nPixDeadSensors)
			self.mu_staco_nSCTDeadSensors.push_back(nSCTDeadSensors)

			self.mu_staco_id_d0.push_back(trackd0)
			self.mu_staco_id_z0.push_back(trackz0)
			self.mu_staco_id_phi.push_back(trackphi)
			self.mu_staco_id_theta.push_back(tracktheta)
			self.mu_staco_id_qoverp.push_back(trackqoverp)

			self.mu_staco_id_theta_exPV.push_back(tracktheta)
			self.mu_staco_id_qoverp_exPV.push_back(trackqoverp)

			self.mu_staco_trackd0pvunbiased.push_back(trackd0pvunbiased)
			self.mu_staco_trackz0pvunbiased.push_back(trackz0pvunbiased)
			self.mu_staco_tracksigd0pvunbiased.push_back(tracksigd0pvunbiased)
			self.mu_staco_tracksigz0pvunbiased.push_back(tracksigz0pvunbiased)

			self.mu_staco_ptcone20.push_back(ptcone20)
			self.mu_staco_ptcone30.push_back(ptcone30)
			self.mu_staco_ptcone40.push_back(ptcone40)

			self.mu_staco_etcone20.push_back(Etcone20)
			self.mu_staco_etcone30.push_back(Etcone30)
			self.mu_staco_etcone40.push_back(Etcone40)

			self.mu_staco_truth_type.push_back(truth_type)
			self.mu_staco_truth_barcode.push_back(truth_barcode)
			self.mu_staco_truth_mothertype.push_back(truth_mothertype)
			self.mu_staco_truth_motherbarcode.push_back(truth_motherbarcode)

			self.mu_staco_EFCB_dr.push_back(EF_dr)
			self.mu_staco_EFCB_index.push_back(EF_index)

			##

			self.mu_staco_n[0] += 1

		#############################################################
		# TRIGGERS						    #
		#############################################################

		if self.TrigDecisionTool.getChainGroup('EF_e60_loose').isPassed():
			self.EF_e60_loose[0] = True
		else:
			self.EF_e60_loose[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_e20_medium').isPassed():
			self.EF_e20_medium[0] = True
		else:
			self.EF_e20_medium[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_e10_medium_mu6').isPassed():
			self.EF_e10_medium_mu6[0] = True
		else:
			self.EF_e10_medium_mu6[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_e10_medium_mu10').isPassed():
			self.EF_e10_medium_mu10[0] = True
		else:
			self.EF_e10_medium_mu10[0] = False

		##

		if self.TrigDecisionTool.getChainGroup('EF_2e10_medium').isPassed():
			self.EF_2e10_medium[0] = True
		else:
			self.EF_2e10_medium[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_2e12_medium').isPassed():
			self.EF_2e12_medium[0] = True
		else:
			self.EF_2e12_medium[0] = False

		##

		if self.TrigDecisionTool.getChainGroup('EF_2g15_loose').isPassed():
			self.EF_2g15_loose[0] = True
		else:
			self.EF_2g15_loose[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_2g20_loose').isPassed():
			self.EF_2g20_loose[0] = True
		else:
			self.EF_2g20_loose[0] = False

		##

		if self.TrigDecisionTool.getChainGroup('EF_mu18').isPassed():
			self.EF_mu18[0] = True
		else:
			self.EF_mu18[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_mu18_MG').isPassed():
			self.EF_mu18_MG[0] = True
		else:
			self.EF_mu18_MG[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_mu20').isPassed():
			self.EF_mu20[0] = True
		else:
			self.EF_mu20[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_mu20_MG').isPassed():
			self.EF_mu20_MG[0] = True
		else:
			self.EF_mu20_MG[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_mu22').isPassed():
			self.EF_mu22[0] = True
		else:
			self.EF_mu22[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_mu22_MG').isPassed():
			self.EF_mu22_MG[0] = True
		else:
			self.EF_mu22_MG[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_mu20i').isPassed():
			self.EF_mu20i[0] = True
		else:
			self.EF_mu20i[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_mu40_MSonly').isPassed():
			self.EF_mu40_MSonly[0] = True
		else:
			self.EF_mu40_MSonly[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_2mu10_loose').isPassed():
			self.EF_2mu10_loose[0] = True
		else:
			self.EF_2mu10_loose[0] = False

		if self.TrigDecisionTool.getChainGroup('EF_mu15_mu10_EFFS').isPassed():
			self.EF_mu15_mu10_EFFS[0] = True
		else:
			self.EF_mu15_mu10_EFFS[0] = False

		#############################################################
		# TRIGGER ELECTRONS					    #
		#############################################################

		electrons = self.StoreGateSvc['HLT_egamma_Electrons']

		#############################################################

		for electron in electrons:
			eta = electron.eta()
			phi = electron.phi()

			self.trig_EF_el_eta.push_back(eta)
			self.trig_EF_el_phi.push_back(phi)

			self.trig_EF_el_EF_e60_loose.push_back(self.isFlagged(eta, phi, 'EF_e60_loose'))
			self.trig_EF_el_EF_e20_medium.push_back(self.isFlagged(eta, phi, 'EF_e20_medium'))
			self.trig_EF_el_EF_e10_medium_mu6.push_back(self.isFlagged(eta, phi, 'EF_e10_medium_mu6'))
			self.trig_EF_el_EF_e10_medium_mu10.push_back(self.isFlagged(eta, phi, 'EF_e10_medium_mu10'))

			self.trig_EF_el_EF_2e10_medium.push_back(self.isFlagged(eta, phi, 'EF_2e10_medium'))
			self.trig_EF_el_EF_2e12_medium.push_back(self.isFlagged(eta, phi, 'EF_2e12_medium'))

			self.trig_EF_el_EF_2g15_loose.push_back(self.isFlagged(eta, phi, 'EF_2g15_loose'))
			self.trig_EF_el_EF_2g20_loose.push_back(self.isFlagged(eta, phi, 'EF_2g20_loose'))

			##

			self.trig_EF_el_n[0] += 1

		#############################################################
		# TRIGGER MUONS						    #
		#############################################################

		muons = self.StoreGateSvc['HLT_MuonEFInfo']

		#############################################################

		for muon in muons:
			eta = muon.ExtrapolatedTrack().eta()
			phi = muon.ExtrapolatedTrack().phi()

			if muon.hasExtrapolatedTrack():
				self.trig_EF_trigmuonef_EF_mu18.push_back(self.isFlagged(eta, phi, 'EF_mu18'))
				self.trig_EF_trigmuonef_EF_mu18_MG.push_back(self.isFlagged(eta, phi, 'EF_mu18_MG'))
				self.trig_EF_trigmuonef_EF_mu20.push_back(self.isFlagged(eta, phi, 'EF_mu20'))
				self.trig_EF_trigmuonef_EF_mu20_MG.push_back(self.isFlagged(eta, phi, 'EF_mu20_MG'))
				self.trig_EF_trigmuonef_EF_mu22.push_back(self.isFlagged(eta, phi, 'EF_mu22'))
				self.trig_EF_trigmuonef_EF_mu22_MG.push_back(self.isFlagged(eta, phi, 'EF_mu22_MG'))

				self.trig_EF_trigmuonef_EF_mu20i.push_back(self.isFlagged(eta, phi, 'EF_mu20i'))

				self.trig_EF_trigmuonef_EF_2mu10_loose.push_back(self.isFlagged(eta, phi, 'EF_2mu10_loose'))
				self.trig_EF_trigmuonef_EF_mu15_mu10_EFFS.push_back(self.isFlagged(eta, phi, 'EF_mu15_mu10_EFFS'))
			else:
				self.trig_EF_trigmuonef_EF_mu18.push_back(False)
				self.trig_EF_trigmuonef_EF_mu18_MG.push_back(False)
				self.trig_EF_trigmuonef_EF_mu20.push_back(False)
				self.trig_EF_trigmuonef_EF_mu20_MG.push_back(False)
				self.trig_EF_trigmuonef_EF_mu22.push_back(False)
				self.trig_EF_trigmuonef_EF_mu22_MG.push_back(False)

				self.trig_EF_trigmuonef_EF_mu20i.push_back(False)

				self.trig_EF_trigmuonef_EF_2mu10_loose.push_back(False)
				self.trig_EF_trigmuonef_EF_mu15_mu10_EFFS.push_back(False)

			##

			self.trig_EF_trigmuonef_n[0] += 1

		#############################################################
		#############################################################

		self.Tree1.Fill()
		if isEGamma:
			self.Tree2.Fill()

		return PyAthena.StatusCode.Success

#############################################################################

	def isFlagged(self, eta, phi, chain):

		return particleMatching(eta, phi, self.TrigDecisionTool.features(chain).get('TrigRoiDescriptor')('initialRoI'))[0] >= 0

#############################################################################

from AthenaCommon.AlgSequence import AlgSequence

job = AlgSequence()
job += uD3PD('uD3PD')

theApp.EvtMax = 1000

#############################################################################

