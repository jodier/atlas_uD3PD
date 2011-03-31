#############################################################################
# USER OPTIONS								    #
#############################################################################

AtlGeo = 'ATLAS-GEO-16-00-00'

CondDB = 'OFLCOND-SDR-BS7T-04-08'

#############################################################################

InputFormat = 'AOD'

isMC = True

#############################################################################

InputFiles = [
	'AOD.280342._000152.pool.root'
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
	athenaCommonFlags.PoolESDInput.set_Value_and_Lock(InputFiles)
if InputFormat == 'AOD':
	athenaCommonFlags.PoolAODInput.set_Value_and_Lock(InputFiles)

#############################################################################

import AthenaPoolCnvSvc.ReadAthenaPool

ServiceMgr.EventSelector.InputCollections = InputFiles

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

theAtlasExtrapolator = AtlasExtrapolator()
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
		PyCintex.loadDict("libuD3PDDict.so")

		self.StoreGateSvc = PyAthena.py_svc('StoreGateSvc')
		self.THistSvc = PyAthena.py_svc('THistSvc')

		self.Tree1 = self.THistSvc['/%s/egamma'        % Stream] = ROOT.TTree('egamma'       , 'egamma'       )
		self.Tree2 = self.THistSvc['/%s/egammaTrigDec' % Stream] = ROOT.TTree('egammaTrigDec', 'egammaTrigDec')

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
		self.Tree2.Branch('RunNumber', self.RunNumber, 'RunNumber/i')

		self.Tree1.Branch('EventNumber', self.EventNumber, 'EventNumber/i')
		self.Tree2.Branch('EventNumber', self.EventNumber, 'EventNumber/i')

		self.Tree1.Branch('lbn', self.lbn, 'lbn/i')
		self.Tree2.Branch('lbn', self.lbn, 'lbn/i')

		#########################
		# PRIMARY VERTICES	#
		#########################

		self.vxp_n = array.array('i', [0])

		##

		self.Tree1.Branch('vxp_n', self.vxp_n, 'vxp_n/I')

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

		self.el_Etcone20_corrected = ROOT.std.vector(float)()
		self.el_Etcone30_corrected = ROOT.std.vector(float)()
		self.el_Etcone40_corrected = ROOT.std.vector(float)()

		self.el_trackd0pvunbiased = ROOT.std.vector(float)()
		self.el_tracksigd0pvunbiased = ROOT.std.vector(float)()

		self.el_truth_type = ROOT.std.vector(int)()
		self.el_truth_mothertype = ROOT.std.vector(int)()
		self.el_truth_barbode = ROOT.std.vector(int)()
		self.el_truth_motherbarcode = ROOT.std.vector(int)()

		self.el_type = ROOT.std.vector(int)()
		self.el_origin = ROOT.std.vector(int)()
		self.el_typebkg = ROOT.std.vector(int)()
		self.el_originbkg = ROOT.std.vector(int)()

		#

		self.Tree1.Branch('el_n', self.el_n, 'el_n/I')

		self.Tree1.Branch('el_E', self.el_E)
		self.Tree1.Branch('el_Et', self.el_Et)
		self.Tree1.Branch('el_pt', self.el_pt)
		self.Tree1.Branch('el_eta', self.el_eta)
		self.Tree1.Branch('el_phi', self.el_phi)
		self.Tree1.Branch('el_charge', self.el_charge)
		self.Tree1.Branch('el_author', self.el_author)
		self.Tree1.Branch('el_isEM', self.el_isEM)

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

		self.Tree1.Branch('el_Etcone20_corrected', self.el_Etcone20_corrected)
		self.Tree1.Branch('el_Etcone30_corrected', self.el_Etcone30_corrected)
		self.Tree1.Branch('el_Etcone40_corrected', self.el_Etcone40_corrected)

		self.Tree1.Branch('el_trackd0pvunbiased', self.el_trackd0pvunbiased)
		self.Tree1.Branch('el_tracksigd0pvunbiased', self.el_tracksigd0pvunbiased)

		self.Tree1.Branch('el_truth_type', self.el_truth_type)
		self.Tree1.Branch('el_truth_mothertype', self.el_truth_mothertype)
		self.Tree1.Branch('el_truth_barbode', self.el_truth_barbode)
		self.Tree1.Branch('el_truth_motherbarcode', self.el_truth_motherbarcode)

		self.Tree1.Branch('el_type', self.el_type)
		self.Tree1.Branch('el_origin', self.el_origin)
		self.Tree1.Branch('el_typebkg', self.el_typebkg)
		self.Tree1.Branch('el_originbkg', self.el_originbkg)

		#########################
		# MUONS			#
		#########################

		self.mu_muid_n = array.array('i', [0])

		# TODO #

		self.mu_muid_trackd0pvunbiased = ROOT.std.vector(float)()
		self.mu_muid_tracksigd0pvunbiased = ROOT.std.vector(float)()

		# TODO #

		##

		self.Tree1.Branch('mu_muid_n', self.mu_muid_n, 'mu_muid_n/I')

		# TODO #

		self.Tree1.Branch('mu_muid_trackd0pvunbiased', self.mu_muid_trackd0pvunbiased)
		self.Tree1.Branch('mu_muid_tracksigd0pvunbiased', self.mu_muid_tracksigd0pvunbiased)

		# TODO #

		#########################
		# TRIGGERS		#
		#########################

		self.L1_EM10 = array.array('i', [0])
		self.L1_EM14 = array.array('i', [0])

		self.EF_e10_medium = array.array('i', [0])
		self.EF_e15_medium = array.array('i', [0])

		##

		self.Tree2.Branch('L1_EM10', self.L1_EM10, 'L1_EM10/I')
		self.Tree2.Branch('L1_EM14', self.L1_EM14, 'L1_EM14/I')

		self.Tree2.Branch('EF_e10_medium', self.EF_e10_medium, 'EF_e10_medium/I')
		self.Tree2.Branch('EF_e15_medium', self.EF_e15_medium, 'EF_e15_medium/I')

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

		self.el_Etcone20_corrected.clear()
		self.el_Etcone30_corrected.clear()
		self.el_Etcone40_corrected.clear()

		self.el_trackd0pvunbiased.clear()
		self.el_tracksigd0pvunbiased.clear()

		self.el_truth_type.clear()
		self.el_truth_barbode.clear()
		self.el_truth_mothertype.clear()
		self.el_truth_motherbarcode.clear()

		self.el_type.clear()
		self.el_origin.clear()
		self.el_typebkg.clear()
		self.el_originbkg.clear()

		#########################
		# MUONS			#
		#########################

		self.mu_muid_n[0] = 0

		# TODO #

		self.mu_muid_trackd0pvunbiased.clear()
		self.mu_muid_tracksigd0pvunbiased.clear()

		# TODO #

		#########################
		# TRIGGERS		#
		#########################

		self.L1_EM10[0] = 0
		self.L1_EM14[0] = 0

		self.EF_e10_medium[0] = 0
		self.EF_e15_medium[0] = 0

	#####################################################################

	def finalize(self):
		return PyAthena.StatusCode.Success

	#####################################################################

	def execute(self):
		self.treeCleaner()

		#############################################################
		# EVENT							    #
		#############################################################

		event = self.StoreGateSvc['MyEvent']

		eventID = event.event_ID()

		#############################################################

		self.RunNumber[0] = eventID.run_number()
		self.EventNumber[0] = eventID.event_number()
		self.lbn[0] = eventID.lumi_block()

		#############################################################
		# PRIMARY VERTICES					    #
		#############################################################

		vertices = self.StoreGateSvc['VxPrimaryCandidate']

		#############################################################

		self.vxp_n[0] = len(vertices)

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
					if (isEM2 & PyAthena.egammaPID.frwdElectronLoose) == 0:
						loose = True
					if (isEM2 & PyAthena.egammaPID.frwdElectronTight) == 0:
						tight = True
				else:
					if (isEM2 & PyAthena.egammaPID.ElectronLoose) == 0:
						loose = True
					if (isEM2 & PyAthena.egammaPID.ElectronMedium) == 0:
						medium = True
					if (isEM2 & PyAthena.egammaPID.ElectronTight) == 0:
						tight = True

				##

				is_ok = False

				for i in xrange(electron.nDetails()):

					detail = electron.detail(i)

					if detail and isinstance(detail, PyAthena.EMShower) and electron.detailName(i) == 'egDetailAOD':

						reta = detail.e237() /\
							    detail.e277()
						weta2 = detail.weta2()

						ptcone20 = detail.ptcone20()
						ptcone30 = detail.ptcone30()
						ptcone40 = detail.ptcone40()

						Etcone20 = detail.etcone20()
						Etcone30 = detail.etcone30()
						Etcone40 = detail.etcone40()

						Etcone20_corrected = -999999.0
						Etcone30_corrected = -999999.0
						Etcone40_corrected = -999999.0

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

					Es2 = cluster.energyBE(2)
					etas2 = cluster.etaBE(2)
					phis2 = cluster.phiBE(2)

				else:
					continue

				##

				track = electron.trackParticle()

				if track:

					r = self.TrackToVertexIPEstimator.estimate(track, vertices[0], True)

					trackd0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getD0(r)
					tracksigd0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getSigma(r)

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

					epyt = -999999			# TODO #
					origin = -999999		# TODO #
					typebkg = -999999		# TODO #
					originbkg = -999999		# TODO #

				##

				self.el_E.push_back(E)
				self.el_Et.push_back(Et)
				self.el_pt.push_back(pt)
				self.el_eta.push_back(eta)
				self.el_phi.push_back(phi)
				self.el_charge.push_back(charge)
				self.el_author.push_back(author)
				self.el_isEM.push_back(isEM1)

				self.el_cl_E.push_back(cl_E)
				self.el_cl_pt.push_back(cl_pt)
				self.el_cl_eta.push_back(cl_eta)
				self.el_cl_phi.push_back(cl_phi)

				self.el_Es2.push_back(Es2)
				self.el_etas2.push_back(etas2)
				self.el_phis2.push_back(phis2)

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

				self.el_Etcone20_corrected.push_back(Etcone20_corrected)
				self.el_Etcone30_corrected.push_back(Etcone30_corrected)
				self.el_Etcone40_corrected.push_back(Etcone40_corrected)

				self.el_trackd0pvunbiased.push_back(trackd0pvunbiased)
				self.el_tracksigd0pvunbiased.push_back(tracksigd0pvunbiased)

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
		# MUONS							    #
		#############################################################

		muons = self.StoreGateSvc['MuidMuonCollection']

		#############################################################

		if len(vertices) > 0:
			for muon in muons:

				# TODO #

				##

				# TODO #

				##

				track = muon.track()

				if track:

					r = self.TrackToVertexIPEstimator.estimate(track, vertices[0], True)

					trackd0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getD0(r)
					tracksigd0pvunbiased = ROOT.Trk.ImpactParametersAndSigma__getSigma(r)

				else:
					continue

				##

				# TODO #

				##

				# TODO #

				self.mu_muid_trackd0pvunbiased.push_back(trackd0pvunbiased)
				self.mu_muid_tracksigd0pvunbiased.push_back(tracksigd0pvunbiased)

				# TODO #

				##

				self.mu_muid_n[0] += 1

		#############################################################
		# TRIGGERS						    #
		#############################################################

		ctp = self.StoreGateSvc['CTP_Decision']

		ctp_items = ctp.getItems()

		#############################################################

		self.L1_EM10[0] = 0
		self.L1_EM14[0] = 0

		for ctp_item in ctp_items:

			if ctp_item == 'L1_EM10':
				self.L1_EM10[0] = 1

			if ctp_item == 'L1_EM14':
				self.L1_EM14[0] = 1

		if self.TrigDecisionTool.getChainGroup('EF_e10_medium').isPassed():
			self.EF_e10_medium[0] = 1
		else:
			self.EF_e10_medium[0] = 0

		if self.TrigDecisionTool.getChainGroup('EF_e15_medium').isPassed():
			self.EF_e15_medium[0] = 1
		else:
			self.EF_e15_medium[0] = 0

		#############################################################
		#############################################################

		self.Tree1.Fill()
		self.Tree2.Fill()

		return PyAthena.StatusCode.Success

#############################################################################

from AthenaCommon.AlgSequence import AlgSequence

job = AlgSequence()
job += uD3PD('uD3PD')

#############################################################################

