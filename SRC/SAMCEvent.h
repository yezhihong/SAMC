#ifndef RAD_PHYSICAL_CONSTANT_H
#include "SAMC.h"
//#include "Your_Cross_Section.h" //The file you store your own cross section models, call it when caluclating cs_Final
/*Call fortran codes{{{*/
extern "C"
{
	//must be like the following
	//change l(r)txfit_ in fortran to real function ...
	//otherwise no result is passed.
	//for left arm
	float ltxfit_(float*,int&);
	float ldelta_(float*,int&);
	float ltheta_(float*,int&);
	float lphi_(float*,int&);
	float ly00_(float*,int&);
	float x_e_q1ex_(float*,int&);
	float y_e_q1ex_(float*,int&);
	float x_e_dent_(float*,int&);
	float y_e_dent_(float*,int&);
	float x_e_dext_(float*,int&);
	float y_e_dext_(float*,int&);
	float x_e_q3en_(float*,int&);
	float y_e_q3en_(float*,int&);
	float x_e_q3ex_(float*,int&);
	float y_e_q3ex_(float*,int&);
	float x_e_fp_(float*,int&);
	float y_e_fp_(float*,int&);
	float p_e_fp_(float*,int&);
	float t_e_fp_(float*,int&);
	void left_init_r_function_();
	void left_getindex_(float*,float*,int*);
	float left_rfunction_(float*,float*,float*,float*);

	//for right arm
	float rtxfit_(float*,int&);
	float rdelta_(float*,int&);
	float rtheta_(float*,int&);
	float rphi_(float*,int&);
	float ry00_(float*,int&);
	float x_h_q1ex_(float*,int&);
	float y_h_q1ex_(float*,int&);
	float x_h_dent_(float*,int&);
	float y_h_dent_(float*,int&);
	float x_h_dext_(float*,int&);
	float y_h_dext_(float*,int&);
	float x_h_q3en_(float*,int&);
	float y_h_q3en_(float*,int&);
	float x_h_q3ex_(float*,int&);
	float y_h_q3ex_(float*,int&);
	float x_h_fp_(float*,int&);
	float y_h_fp_(float*,int&);
	float p_h_fp_(float*,int&);
	float t_h_fp_(float*,int&);
	void right_init_r_function_();
	void right_getindex_(float*,float*,int*);
	float right_rfunction_(float*,float*,float*,float*);
}
/*}}}*/
#define MSIZE 5
/*struct Material{{{*/
struct Material
{
	string Name; //Name
	int Z;//Z
	double A;//A
	double M;//Mass
	double T;//Thickness g/cm^2
	double TR;//Thickness in Rad_Len
	double rho;//density g/cm^3
	double bt;//bt
	double X0; //Radiation Length g/cm^2
	double L;//Length:distance along Z axix in TCS
}
;/*}}}*/
#endif

/*class SAMCEvent{{{*/
class SAMCEvent
{
	public:
		/*SAMCEvent(){{{*/
		SAMCEvent()
		{
			Init();
		}
		/*}}}*/
		/*virtual ~SAMCEvent(){{{*/
		virtual ~SAMCEvent()
		{
			Win_Before_Mag.clear();
			Win_After_Mag.clear();
		}
		/*}}}*/

		SAMCEvent(SAMCEvent const&){};
		SAMCEvent& operator=(SAMCEvent const&){};

		/*int Process(){{{*/
		int Process()
		{
			int Num_Event_Add=0;
			Generator();//generate original target variables
			Num_Event_Add=RefineTg();//transport to front of magnetic, then back to get refined target variables
			if ( Num_Event_Add>0 ) {
				IsPassed=0;
				IsQualified=0;
			}
			else
				Num_Event_Add=ToFp(x_tg_ref,y_tg_ref,th_tg_ref,ph_tg_ref,dp_ref);//transfer to focal plane using John.LeRose matrix
			if ( Num_Event_Add>0 ) {
				IsPassed=0;
				IsQualified=0;
			}
			else
				Num_Event_Add=ReconstructTg(x_fp,y_fp,th_fp,ph_fp,x_tg_ref);
			return 0;
		}
		/*}}}*/

		/*void AddOneMaterial(vector<Material>& aWin,const double& aX0,const double& arho,const double& aL,const double& aA,const int& aZ,string aName){{{*/
		void AddOneMaterial(vector<Material>& aWin,const double& aX0,const double& arho,const double& aL,const double& aA,const int& aZ,string aName)
		{
			//To make sure read correctly, so I reverse the order compared with input file.
			Material a;
			a.Name=aName;
			a.Z=aZ;
			a.A=aA;
			a.L=aL;
			a.rho=arho;
			a.T=aL*arho;
			a.X0=aX0;
			if ( fabs(aX0)<1e-10 )
			{
				a.TR=0;
				a.bt=0;
			}
			else
			{
				a.TR=a.T/aX0;
				a.bt=b(aZ)*aL*arho/aX0;
			}
			aWin.push_back(a);
		}
		/*}}}*/

	private:
		/*void Generator(){{{*/
		void Generator()
		{
			//set value for Member Data derived from variables from file
			//File provides E_s,theta,Target.(Z,A,T,rho) Win_i.(Z,A,T,rho)
			//Win_f.(Z,A,T,rho) T_theta
			//Win_Before_Mag(Name,Z,A,L,rho,X0)
			//Win_After_Mag(Name,Z,A,L,rho,X0)
			//beam_x,beam_y,reactz_gen,th_tg_gen,ph_tg_gen,dp_gen
			//z0,HRS_L,VDC_Res(x,y,th,ph),D_(x,y),T_L,P0
			//IsMultiScat,IsEnergyLoss,Which_Kin,FP_Eff_L

			/*Set Material.(Z,A,M,X0,T,TR,bt){{{*/
			SetMaterial(Target);
			SetMaterial(Win_i);
			SetMaterial(Win_f);
			/*}}}*/

			/*Set Beam Info(HCS and TCS){{{*/
			//know beam_x,beam_y,reactz_gen,E_s,theta,HRS_L
			//Set s,s_TCS,x_tg_gen,y_tg_gen,p_TCS,p_P,p_P_TCS
			theta_rad=theta*DegToRad();//rad
			s(0)=beam_x; //cm
			s(1)=beam_y; //cm
			reactz_gen+=-(beam_x)*tan(T_theta*DegToRad());
			s(2)=reactz_gen; //cm
			target_edgepoint_TRCS(0)=theta/fabs(theta)*T_H/2;//if theta>0,T_H/2, if<0, -T_H/2 in Target Rotation Coordinate System(T_theta=0) not TCS, check Coordinate.svg
			target_edgepoint_TRCS(1)=0;
			target_edgepoint_TRCS(2)=T_L/2;//and z0=0

			TLorentzVector lp_TRCS;//the interaction point in TRCS at T_theta=0
			lp_TRCS=s;
			lp_TRCS(2)-=z0;
			//Printf("s(%g,%g,%g),target_edgepoint_TRCS(%g,%g,%g),lp_TRCS(%g,%g,%g)",s(0),s(1),s(2),target_edgepoint_TRCS(0),target_edgepoint_TRCS(1),target_edgepoint_TRCS(2),lp_TRCS(0),lp_TRCS(1),lp_TRCS(2));
			lp_TRCS.RotateY(-T_theta*DegToRad());

			s_TCS=s;
			s_TCS.RotateZ(PI/2);//passive ratation around HCS, so -(-PI/2)
			s_TCS.RotateX(theta_rad);//passive ratation around HCS, so -(-theta_rad)
			s_TCS(0)-=D_x;
			s_TCS(1)-=D_y;

			p_inter_point_TCS=s_TCS;
			reactz_TCS=s_TCS.Z();

			s_TCS(0)-=s_TCS.Z()*th_tg_gen;
			s_TCS(1)-=s_TCS.Z()*ph_tg_gen;
			s_TCS(2)-=s_TCS.Z();

			x_tg_gen=s_TCS(0);
			y_tg_gen=s_TCS(1);
			
			//Fix me: I don't think we need this
			//Material halftarget=Target;
			//halftarget.T /= 2;
			//E_s-=Ion_Loss(E_s,Win_i);
			//E_s-=Bremss_Loss(E_s,Win_i.bt+btr);
			//E_s-=Ion_Loss(E_s,halftarget);
			//E_s-=Bremss_Loss(E_s,halftarget.bt+btr);
			s(3)=E_s;//MeV
			s_TCS(3)=E_s;//MeV

			p_TCS.SetVect(s_TCS.Vect());
			TLorentzVector lz(0,0,1,0);//in HCS
			p_P=lz;
			p_P.RotateZ(PI/2);
			p_P.RotateX(theta_rad);
			p_P.RotateY(-atan(th_tg_gen));
			p_P.RotateX(atan(ph_tg_gen));
			Angle=p_P.Angle(lz.Vect());
			Angle_Deg=Angle*RadToDeg();
			sinsq_Angle=sin(Angle/2)*sin(Angle/2);
			sinsq=sin(Angle/2)*sin(Angle/2);//sin(Angle/2)^2
			//p_P_TCS=lz;//Now think it's in TCS
			//p_P_TCS.RotateY(atan(th_tg_gen));//p_P_TCS.X()/p_P_TCS.Z()=th_tg
			//p_P_TCS.RotateX(-atan(ph_tg_gen));//p_P_TCS.Y()/p_P_TCS.Z()=ph_tg
			p_P_TCS.SetX(th_tg_gen);
			p_P_TCS.SetY(ph_tg_gen);
			p_P_TCS.SetZ(1);
			p_P_TCS.SetVect(p_P_TCS.Vect().Unit());
			/*}}}*/

			size_t i;
			size_t imax;
			Material m0;
			/*Set Windows{{{*/
			//know Win_Before_Mag,Win_f,HRS_L,theta
			//Add Target,Win_f,Air if necessary to Win_Before_Mag
			//Win_After_Mag don't need to be changed since it's vacuum.
			//for Ion_Loss Z,A,T,rho
			//for Bremss_Loss bt
			//for MultiScattering TR
			//for Transport L
			//Thickness definition diagram
			//Win_i |<--T-->|

			//Target|<--T-->|

			//Win_f ---------
			//      |<--T
			//      ---------

			//Target diagram,z0,TCS,reactz_gen,T
			//|<----T--->|
			//|  !  !   !|
			//   TCS!   !
			//      z0  !
			//          reactz_gen
			//          <-lL
			vector<Material>::iterator it=Win_Before_Mag.begin();//here only have the material before mag defined in input file.
			double lHL=0; //if lHL<HRS_Lcm or <3.57m(See Hall_A_NIM), i assume it's air
			int FirstInWinBeforeMagBlock=0;
			imax=Win_Before_Mag.size();
			
			for ( i=0; i<imax; i++ )
				lHL+=Win_Before_Mag[i].L;

			double lphrad=theta_rad-T_theta*DegToRad();//scattering angle in y axis in TCS or x axis in TRCS
			//it doesn't I add atan(ph_tg_gen) since it's so small for target and windows



			//Win_Before_Mag add Win_f
			if ( Win_f.Z!=0 && Win_f.A!=0 && Win_f.rho!=0 && Win_f.TR!=0 )
			{
				m0=Win_f;//Win_f added to Win_Before_Mag
				m0.Name="Win_f";
				Win_Before_Mag.insert(it,m0);
				FirstInWinBeforeMagBlock++;
			}
			//Win_Before_Mag add Target
			it=Win_Before_Mag.begin();

			m0=Target;//target after interaction point added to Win_Before_Mag
			double lL; //assume target is rentangle lL=distance between reactz_gen and edge of target
			//in TRCS x plane<--> y plane in TCS
			lL=target_edgepoint_TRCS(2)-lp_TRCS(2);//verticle distance between interaction point and edge plane of target
			if ( lL>T_L ) {
				Printf("target_edgepoint_TRCS(%g,%g,%g),lp_TRCS(%g,%g,%g)",target_edgepoint_TRCS(0),target_edgepoint_TRCS(1),target_edgepoint_TRCS(2),lp_TRCS(0),lp_TRCS(1),lp_TRCS(2));
			}
			lL*=tan(lphrad);
			lL+=lp_TRCS(0);//x on edge
			
			if ( fabs(lL)<fabs(target_edgepoint_TRCS(0)) ) {//HCS 0=x top view
				//in the target
				lL=fabs((target_edgepoint_TRCS(2)-lp_TRCS(2))/cos(lphrad));
			}
			else {
				//out of target before hiting the edge of target, the edge means the downstream face.
				lL=fabs((target_edgepoint_TRCS(0)-lp_TRCS(0))/sin(lphrad));
			}
			
			m0.L=lL;
			Win_Before_Mag.insert(it,m0);
			FirstInWinBeforeMagBlock++;

			//Win_Before_Mag add Air before last material in the input file
			imax=Win_Before_Mag.size();
			if ( (lHL)<(HRS_L) )
			{
				//http://pdg.lbl.gov/2009/AtomicNuclearProperties/HTML_PAGES/104.html
				double total=0.000124+0.755267+0.231781+0.012827;//total mass of component of air
				m0.Name="Air";
				m0.Z=int(6*0.000124/total+7*0.755267/total+8*0.231781/total+18*0.012827/total);
				m0.A=m0.Z/0.499;
				m0.L=HRS_L-lHL;
				
				m0.rho=1.205e-03;
				m0.T=m0.L*m0.rho;
				m0.X0=36.66;
				m0.bt=b(m0.Z)*m0.T/m0.X0;
				it=Win_Before_Mag.end();
				it--;
				Win_Before_Mag.insert(it,m0);
			}
			else if ( (lHL)>HRS_L )
			{
				printf("[Error %s: Line %d] Total windows length=%f>HRS_L=%f.\n",__FILE__,__LINE__,(lHL+Win_Before_Mag[imax-1].L),HRS_L);
				exit(-3);
			}
			//Correct First Material.L in Win_Before_Mag Block of inputfile
			it=Win_Before_Mag.begin();
			(it+FirstInWinBeforeMagBlock)->L+=reactz_TCS;
			for ( i=0; i<FirstInWinBeforeMagBlock; i++ )
				(it+FirstInWinBeforeMagBlock)->L-=(it+i)->L;

			for ( i = 0; i < Win_Before_Mag.size(); ++i ) {
				Win_Before_Mag[i].T=Win_Before_Mag[i].L*Win_Before_Mag[i].rho;
				if ( fabs(Win_Before_Mag[i].X0)<1e-10 ) {
					Win_Before_Mag[i].TR=0;
				}
				else {
					Win_Before_Mag[i].TR=Win_Before_Mag[i].T/Win_Before_Mag[i].X0;
				}
				Win_Before_Mag[i].bt=b(Win_Before_Mag[i].Z)*Win_Before_Mag[i].TR;
			}
			/*}}}*/

			/*Set Target Info(TCS){{{*/
			//know P0,dp_gen,Which_Kin,IsEnergyLoss,E_s,E_p,sinsq,Target
			//Set dp_gen,s_TCS,x_tg_gen,y_tg_gen,p_TCS,p_P,p_P_TCS
			//Set Q2,q2,btr,Win_Before_Mag[0].bt
			//x_tg_gen,y_tg_gen,see Set Beam Info(HCS and TCS)
			switch ( Which_Kin )
			{
				case 1: //elastic
					E_p=gRandom->Gaus(0,P0*3e-4);
					E_p+=E_s/(1+2*E_s*sinsq/Target.M);
					break;
				case 2: //quasi-elastic
				case 0: //phase default
				default:
					E_p=P0*(1+dp_gen);
			}
			Q2=4*E_s*E_p*sinsq;
			q2=-Q2;
			btr=AP*(log(Q2/(ELECTRON_MASS*ELECTRON_MASS))-1);//b*t_r
			////Correct Win_Before_Mag[0](Target) bt
			//Win_Before_Mag[0].bt+=btr;

			p_TCS(3)=E_p;//MeV
			p_P(3)=E_p;
			p_P_TCS(3)=E_p;
			dp_gen=(E_p-P0)/P0;
			/*}}}*/

		}
		/*}}}*/

		/*int RefineTg(){{{*/
		int RefineTg()
		{
			//return Number of Event needs to be added. just 1 ^_^
			//refine target variables from Generator because John.LeRose matrix only works for vacuum
			//know Win_Before_Mag,p_TCS,p_P_TCS
			p_TCS_ref=p_inter_point_TCS;
			p_P_TCS_ref=p_P_TCS;

			size_t i;
			size_t imax;
			vector<Material> Win_Empty;//tmp use

			Material mixture;

			double offset=p_TCS(2)-p_TCS_ref(2);//L=along Z in TCS
			mixture.L=offset;

			//Printf("p_TCS(%g,%g,%g),(th=%g,ph=%g)",p_TCS(0),p_TCS(1),p_TCS(2),th_tg_gen,ph_tg_gen);
			//FIXME: I want to add energy loss in GetRef_Plane
			//but the SAMC will run forever and no good event.
			if ( E_p<ELECTRON_MASS ) {
				printf("[Warning %s: Line %d] E_p=%g<ELECTRON_MASS=%g\n",__FILE__,__LINE__,E_p,ELECTRON_MASS);
				return 1;
			}

			GetRef_Plane(p_TCS_ref,p_P_TCS_ref,Win_Before_Mag,Win_Empty,0,offset);

			x_tg_ref=p_TCS_ref(0);
			y_tg_ref=p_TCS_ref(1);
			th_tg_ref=p_P_TCS_ref(0)/p_P_TCS_ref(2);
			ph_tg_ref=p_P_TCS_ref(1)/p_P_TCS_ref(2);
			//Printf("p_TCS_ref(%g,%g,%g),(th=%g,ph=%g)",p_TCS_ref(0),p_TCS_ref(1),p_TCS_ref(2),th_tg_ref,ph_tg_ref);

			dp_ref=dp_gen;
			E_p=p_P_TCS_ref(3);

			if ( IsEnergyLoss )
			{
				//Fix me: Do I need to count the energy loss due to windows after target?
				//I think so. So I add them.
				imax=Win_Before_Mag.size();
				for ( i=0; i<imax; i++ )
				{
					E_p-=Ion_Loss(E_p,Win_Before_Mag[i]);
					E_p-=Bremss_Loss(E_p,Win_Before_Mag[i].bt);
				}
				imax=Win_After_Mag.size();
				for ( i=0; i<imax; i++ )
				{
					E_p-=Ion_Loss(E_p,Win_After_Mag[i]);
					E_p-=Bremss_Loss(E_p,Win_After_Mag[i].bt);
				}
				if ( E_p<0 ) {
					printf("[Warning %s: Line %d] E_p after Energy Loss=%f<0.\n",__FILE__,__LINE__,E_p);
					return 1;
				}
				Q2=4*E_s*E_p*sinsq;
				q2=-Q2;
				p_TCS_ref(3)=E_p;//MeV
				p_P_TCS_ref(3)=E_p;
				dp_ref=(E_p-P0)/P0;
			}
			
			
			cs_M=sigma_M(E_s,Angle_Deg);
			cs_Final=cs_M;//Just use the Mott cross section so far
			//////////////////////////////////////////////////////////////////
			// This is where you add your own cross section.
			// e.g.:
			// cs_Final = Your_Cross_Section(E_s, E_p, Angle_Deg, A, Z);
            //////////////////////////////////////////////////////////////////
			return 0;
		}
		/*}}}*/

		/*int ToFp(const double& ax,const double& ay,const double& ath,const double& aph,const double& adp){{{*/
		int ToFp(const double& ax,const double& ay,const double& ath,const double& aph,const double& adp)
		{
			//return Number of Event needs to be added. just 1 ^_^
			//cm->meter, ath,aph,adp are correct
			//must be float, otherwise cannot pass to fortran correctly. float<->dimension
			float matrix[MSIZE]={ax/100.,ath,ay/100,aph,adp};
			int msize=MSIZE;
			double xtest,ytest;
			unsigned int i;
			for ( i = 0; i < 2; ++i ) {
				q1ex[i]=0;
				dent[i]=0;
				dext[i]=0;
				q3en[i]=0;
				q3ex[i]=0;
			}
			IsPassedQ1Ex=false;
			IsPassedDipoleEn=false;
			IsPassedDipoleEx=false;
			IsPassedQ3En=false;
			IsPassedQ3Ex=false;
			if ( theta>0 )
			{
				//for left arm
				/*Q1 Exit{{{*/
				xtest=x_e_q1ex_(matrix,msize)*100; //cm
				ytest=y_e_q1ex_(matrix,msize)*100; //cm
				q1ex[0]=xtest/100;
				q1ex[1]=ytest/100;
				if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q1_Radius*Q1_Radius )
				{
					//(xtets!=xtest)==true to avoid nan(Not a number)
					//printf("Blocked by Q1 Exit.\n");
					return 1;
				}
				IsPassedQ1Ex=true;
				/*}}}*/
				/*Dipole Entrance{{{*/
				// Transport electron to dipole entrance, trapezoid define by (jjl)
				// -40cm<x<40cm (+x is down in HRS frame)
				// y=+-(12.5*(1-(1.25*x/840)) (smallest gap at positive x, x in cm)
				xtest=x_e_dent_(matrix,msize)*100; //cm
				ytest=y_e_dent_(matrix,msize)*100; //cm
				dent[0]=xtest/100;
				dent[1]=ytest/100;
				if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
				{
					//printf("Blocked by Dipole Entrance.\n");
					return 1;
				}
				IsPassedDipoleEn=true;
				/*}}}*/

				/*Dipole Exit{{{*/
				xtest=x_e_dext_(matrix,msize)*100; //cm
				ytest=y_e_dext_(matrix,msize)*100; //cm
				dext[0]=xtest/100;
				dext[1]=ytest/100;
				if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
				{
					//printf("Blocked by Dipole Exit.\n");
					return 1;
				}
				IsPassedDipoleEx=true;
				/*}}}*/

				/*Q3 Entrance{{{*/
				xtest=x_e_q3en_(matrix,msize)*100; //cm
				ytest=y_e_q3en_(matrix,msize)*100; //cm
				q3en[0]=xtest/100;
				q3en[1]=ytest/100;
				if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q3_Entrance_Radius*Q3_Entrance_Radius )
				{
					//printf("Blocked by Q3 Entrance.\n");
					return 1;
				}
				IsPassedQ3En=true;

				/*}}}*/
 
				/*Q3 Exit{{{*/
				xtest=x_e_q3ex_(matrix,msize)*100; //cm
				ytest=y_e_q3ex_(matrix,msize)*100; //cm
				q3ex[0]=xtest/100;
				q3ex[1]=ytest/100;
				if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q3_Exit_Radius*Q3_Exit_Radius )
				{
					//printf("Blocked by Q3 Exit.\n");
					return 1;
				}
				IsPassedQ3Ex=true;
				/*}}}*/

				x_fp=x_e_fp_(matrix,msize)*100.; //cm
				y_fp=y_e_fp_(matrix,msize)*100.; //cm
				th_fp=t_e_fp_(matrix,msize); //(tantheta)
				ph_fp=p_e_fp_(matrix,msize); //(tanphi)
			}
			else
			{
				//for right arm
				/*Q1 Exit{{{*/
				xtest=x_h_q1ex_(matrix,msize)*100; //cm
				ytest=y_h_q1ex_(matrix,msize)*100; //cm
				q1ex[0]=xtest/100;
				q1ex[1]=ytest/100;
				if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q1_Radius*Q1_Radius )
				{
					//(xtets!=xtest)==true to avoid nan(Not a number)
					//printf("Blocked by Q1 Exit.\n");
					return 1;
				}
				IsPassedQ1Ex=true;
				/*}}}*/
	
				/*Dipole Entrance{{{*/
				// Transport electron to dipole entrance, trapezoid define by (jjl)
				// -40cm<x<40cm (+x is down in HRS frame)
				// y=+-(D_Y_L*(1-(1.25*x/840)) (smallest gap at positive x, x in cm)
				xtest=x_h_dent_(matrix,msize)*100; //cm
				ytest=y_h_dent_(matrix,msize)*100; //cm
				dent[0]=xtest/100;
				dent[1]=ytest/100;
				if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
				{
					//printf("Blocked by Dipole Entrance.\n");
					return 1;
				}
				IsPassedDipoleEn=true;
	 			/*}}}*/
	
				/*Dipole Exit{{{*/
				xtest=x_h_dext_(matrix,msize)*100; //cm
				ytest=y_h_dext_(matrix,msize)*100; //cm
				dext[0]=xtest/100;
				dext[1]=ytest/100;
				if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
				{
					//printf("Blocked by Dipole Exit.\n");
					return 1;
				}
				IsPassedDipoleEx=true;
				/*}}}*/
	
				/*Q3 Entrance{{{*/
				xtest=x_h_q3en_(matrix,msize)*100; //cm
				ytest=y_h_q3en_(matrix,msize)*100; //cm
				q3en[0]=xtest/100;
				q3en[1]=ytest/100;
				if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q3_Entrance_Radius*Q3_Entrance_Radius )
				{
					//printf("Blocked by Q3 Entrance.\n");
					return 1;
				}
				IsPassedQ3En=true;
				/*}}}*/
		
				/*Q3 Exit{{{*/
				xtest=x_h_q3ex_(matrix,msize)*100; //cm
				ytest=y_h_q3ex_(matrix,msize)*100; //cm
				q3ex[0]=xtest/100;
				q3ex[1]=ytest/100;
				if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q3_Exit_Radius*Q3_Exit_Radius )
				{
					//printf("Blocked by Q3 Exit.\n");
					return 1;
				}
				IsPassedQ3Ex=true;
				/*}}}*/

				x_fp=x_h_fp_(matrix,msize)*100.; //cm
				y_fp=y_h_fp_(matrix,msize)*100.; //cm
				th_fp=t_h_fp_(matrix,msize); //(tantheta)
				ph_fp=p_h_fp_(matrix,msize); //(tanphi)
			}
			if ( IsMultiScat)
			{
				//Material mixture=GetMixture(Win_After_Mag);
				//mixture.TR*=sqrt(ph_fp*ph_fp+th_fp*th_fp);
				//mixture.L=0;
				//p_FP.SetXYZT(x_fp,y_fp,0,p_TCS_ref(3));
				//p_P_FP.SetXYZT(th_fp,ph_fp,1,p_TCS_ref(3));
				//Transport(p_FP,p_P_FP,mixture,IsMultiScat);//Just change th_fp,ph_fp
				//th_fp=p_P_FP(0)/p_P_FP(2);
				//ph_fp=p_P_FP(1)/p_P_FP(2);

				//Material mixture=GetMixture(Win_After_Mag);
				p_FP.SetXYZT(x_fp,y_fp,0,p_TCS_ref(3));
				p_P_FP.SetXYZT(th_fp,ph_fp,1,p_TCS_ref(3));
				Material Vacuum;
				Vacuum.L=0;
				for ( i = 0; i < Win_After_Mag.size(); ++i ) {
					Vacuum.L-=Win_After_Mag[i].L;
				}
				Vacuum.Name="Vacuum";
				Vacuum.TR=0;
				Transport(p_FP,p_P_FP,Vacuum,IsMultiScat);//Just change th_fp,ph_fp
				for ( i = 0; i < Win_After_Mag.size(); ++i ) {
					Transport(p_FP,p_P_FP,Win_After_Mag[i],IsMultiScat);//Just change th_fp,ph_fp
				}
				//Transport(p_FP,p_P_FP,mixture,IsMultiScat);//Just change th_fp,ph_fp
				th_fp=p_P_FP(0)/p_P_FP(2);
				ph_fp=p_P_FP(1)/p_P_FP(2);
				x_fp=p_FP(0);
				y_fp=p_FP(1);
			}
			//Smearing
			TRandom* tr=new TRandom();
			x_fp=tr->Gaus(x_fp,VDC_Res_x);
			y_fp=tr->Gaus(y_fp,VDC_Res_y);
			th_fp=tr->Gaus(th_fp,VDC_Res_th/1000.);
			ph_fp=tr->Gaus(ph_fp,VDC_Res_ph/1000.);
			th_fp_no_ort=th_fp;
			delete tr;
			matrix[0]=x_fp/100.;
			msize=1;
			//orthogonalize theta see John.LeRose Webpage
			//http://hallaweb.jlab.org/news/minutes/tranferfuncs.html
			if ( theta>0 )
			{
				th_fp-=ltxfit_(matrix,msize);
			}
			else
				th_fp-=rtxfit_(matrix,msize);
			return 0;
		}
		/*}}}*/

		/*int ReconstructTg(const double& ax,const double& ay,const double& ath,const double& aph,const double& axtg){{{*/
		int ReconstructTg(const double& ax,const double& ay,const double& ath,const double& aph,const double& axtg)
		{
			//return Number of Event needs to be added. just 1 ^_^
			float matrix[MSIZE]={ax/100.,ath,ay/100,aph,axtg/100};
			int msize=MSIZE;

			x_tg_rec=axtg;
			float rf_y,rf_d,rf_th,rf_ph;
			if ( theta>0 )
			{
				//printf("%f %f %f %f %f\n",matrix[0],matrix[1],matrix[2],matrix[3],matrix[4]);
				y_tg_rec=ly00_(matrix,msize)*100;
				th_tg_rec=ltheta_(matrix,msize);
				ph_tg_rec=lphi_(matrix,msize);
				dp_rec=ldelta_(matrix,msize);
				rf_y=y_tg_rec/100;
				rf_d=dp_rec;
				rf_th=th_tg_rec;
				rf_ph=ph_tg_rec;
				if ( fNCuts<=0 ) {
					rvalue=left_rfunction_(&rf_y,&rf_d,&rf_th,&rf_ph);
				}
			}
			else
			{
				y_tg_rec=ry00_(matrix,msize)*100;
				th_tg_rec=rtheta_(matrix,msize);
				ph_tg_rec=rphi_(matrix,msize);
				dp_rec=rdelta_(matrix,msize);
				rf_y=y_tg_rec/100;
				rf_d=dp_rec;
				rf_th=th_tg_rec;
				rf_ph=ph_tg_rec;
				if ( fNCuts<=0 ) {
					rvalue=right_rfunction_(&rf_y,&rf_d,&rf_th,&rf_ph);
				}
			}
			if ( fNCuts>0 ) {
				rvalue=CalcRValue(rf_th,rf_ph,rf_y,rf_d);
			}
			TLorentzVector rec(x_tg_rec,y_tg_rec,0,P0*(1+dp_rec));
			
			rec(0)+=D_x;
			rec(1)+=D_y;
			rec(0)+=reactz_TCS*th_tg_rec;
			rec(1)+=reactz_TCS*ph_tg_rec;
			rec(2)+=reactz_TCS;
			rec.RotateX(-theta_rad);
			rec.RotateZ(-PI/2);
			reactz_rec=rec(2);

			Angle_rec = acos( (cos(theta_rad)-ph_tg_rec*sin(theta_rad)) 
					        / sqrt(1.0+pow(th_tg_rec,2)+pow(ph_tg_rec,2)) );
			// cerr <<"New Angle is " << Angle_rec*RadToDeg() <<endl; 
			Qsq = 4.0*E_s*E_p*sin(Angle_rec/2.0)*sin(Angle_rec/2.0);
			//	cerr <<"Qsq is " << Qsq <<endl;
			Xbj = Qsq/2.0/(E_s-E_p)/PROTON_MASS;
			//	cerr <<"Xbj is " << Xbj <<endl;

			if (
					//fabs(reactz_rec/100)<0.1 &&
					//fabs(y_tg_rec/100)<=0.01 &&
					fabs(dp_rec)<delta_dp/2 &&
					fabs(th_tg_rec)<delta_th/2 &&
					fabs(ph_tg_rec)<delta_ph/2
					//fabs(ph_tg_rec)<delta_ph/2
			   )
				IsQualified=1;
			else
				IsQualified=0;
			return 0;
		}
		/*}}}*/

		/*Bool_t IntersectPlaneWithRay( const TVector3& xax,{{{*/
		Bool_t IntersectPlaneWithRay( const TVector3& xax,
				const TVector3& yax,
				const TVector3& org,
				const TVector3& ray_start,
				const TVector3& ray_vect,
				Double_t& length,
				TVector3& intersect )
		{
			// Find intersection point of plane (given by 'xax', 'yax', 'org') with
			// ray (given by 'ray_start', 'ray_vect'). 
			// Returns true if intersection found, else false (ray parallel to plane).
			// Output is in 'length' and 'intersect', where
			//   intersect = ray_start + length*ray_vect
			// 'length' and 'intersect' must be provided by the caller.

			// Calculate explicitly for speed.

			Double_t nom[9], den[9];
			nom[0] = den[0] = xax.X();
			nom[3] = den[3] = xax.Y();
			nom[6] = den[6] = xax.Z();
			nom[1] = den[1] = yax.X();
			nom[4] = den[4] = yax.Y();
			nom[7] = den[7] = yax.Z();
			den[2] = -ray_vect.X();
			den[5] = -ray_vect.Y();
			den[8] = -ray_vect.Z();

			Double_t det1 = den[0]*(den[4]*den[8]-den[7]*den[5])
				-den[3]*(den[1]*den[8]-den[7]*den[2])
				+den[6]*(den[1]*den[5]-den[4]*den[2]);
			if( fabs(det1) < 1e-5 )
				return false;

			nom[2] = ray_start.X()-org.X();
			nom[5] = ray_start.Y()-org.Y();
			nom[8] = ray_start.Z()-org.Z();
			Double_t det2 = nom[0]*(nom[4]*nom[8]-nom[7]*nom[5])
				-nom[3]*(nom[1]*nom[8]-nom[7]*nom[2])
				+nom[6]*(nom[1]*nom[5]-nom[4]*nom[2]);

			length = det2/det1;
			intersect = ray_start + length*ray_vect;
			return true;
		}
		/*}}}*/

	private:
		/*double Ion_Loss(const double& aE0,const Material& aMaterial){{{*/
		double Ion_Loss(const double& aE0,const Material& aMaterial)
		{
			//aT: g/cm^2, arho: g/cm^3
			//Particle Booklet Equ(27.9)
			//printf("Z=%d,A=%f,T=%f,rho=%f\n",aMaterial.Z,aMaterial.A,aMaterial.T,aMaterial.rho);
			double lK=0.307075;// cm^2/g for A=1 g/mol
			double lbetasq=1-ELECTRON_MASS*ELECTRON_MASS/(aE0*aE0);
			double lxi=lK/2*aMaterial.Z/aMaterial.A*aMaterial.T/lbetasq;//aT: g/cm^2
			double lhbarwsq=28.816*28.816*aMaterial.rho*aMaterial.Z/aMaterial.A*1e-12;//MeV arho is density of absorber
			double j=0.200;
			double Delta_p=lxi*(log(2*ELECTRON_MASS*lxi/lhbarwsq)+j);
			double lw=4*lxi;
			double result=0;
			if ( aMaterial.Z!=0 && aMaterial.A!=0 && aMaterial.T!=0 && aMaterial.rho!=0 )
				result=gRandom->Landau(Delta_p,lw);
			if ( result>(aE0-ELECTRON_MASS) )
				result=aE0-ELECTRON_MASS;
			if ( result<0 )
				result=0;
			return result;
		}
		/*}}}*/

		/*double Bremss_Loss(const double& aE0,const double& abt){{{*/
		double Bremss_Loss(const double& aE0,const double& abt)
		{
			//Bremsstrahlung Energy Loss for external and internal(equivalent radiator)
			//Xiaodong Jiang, PhD.thesis Equ 5.15
			//http://filburt.mit.edu/oops/Html/Pub/theses/xjiang.ps
			//*0.999 to avoid lose all energy
			double result=0;
			if ( abt!=0 )
				result=aE0*pow(gRandom->Rndm()*0.999,1./abt);
			if ( result>(aE0-ELECTRON_MASS) )
				result=aE0-ELECTRON_MASS;
			if ( result<0 )
				result=0;
			return result;
		}
		/*}}}*/
		
		/*double eta(const int& aZ){{{*/
		double eta(const int& aZ)
		{
			//Phys.Rev.D 12,1884 A46
			return log(1440*pow(aZ,-2/3.))/log(183*pow(aZ,-1/3.));
		}
		/*}}}*/

		/*double b(const int& aZ){{{*/
		double b(const int& aZ)
		{
			//Phys.Rev.D 12,1884 A45
			if ( aZ!=0 )
				return 4./3.*( 1+1./9.*(aZ+1)/(aZ+eta(aZ))/log(183*pow(aZ,-1/3.)) );
			else
				return 0;
		}
		/*}}}*/

		/*double Rad_Len(const int& aZ,const double& aA){{{*/
		double Rad_Len(const int& aZ,const double& aA)
		{
			//particle book equation 27.20
			//Lrad=elastic form factor F_el, scattering on the nucleus
			//Lrad_prime=inelastic form factor F_inel, scattering on the shell electrons
			//f(Z)=Coulomb correction
			double Lrad,Lrad_prime,f_Z;
			if ( aZ==1 )
			{
				Lrad=5.31;
				Lrad_prime=6.144;
			}
			else if ( aZ==2 )
			{
				Lrad=4.79;
				Lrad_prime=5.621;
			}
			else if ( aZ==3 )
			{
				Lrad=4.74;
				Lrad_prime=5.805;
			}
			else if ( aZ==4 )
			{
				Lrad=4.71;
				Lrad_prime=5.924;
			}
			else
			{
				Lrad=log(184.15*pow(aZ,-1./3.));
				Lrad_prime=log(1194.*pow(aZ,-2./3.));
			}
			double a=ALPHA*aZ;
			a=a*a;
			f_Z=a*(1./(1+a)+0.20206-0.0369*a+0.0083*a*a-0.002*a*a*a);
			if ( aZ!=0 && aA!=0 )
				return 716.408*aA/(aZ*aZ*(Lrad-f_Z)+aZ*Lrad_prime);
			else
				return 0;
		}
		/*}}}*/

		/*void Transport(TLorentzVector& aPos,TLorentzVector& aMom,const Material& aMaterial,const bool& aIsMultiScatt=false){{{*/
		void Transport(TLorentzVector& aPos,TLorentzVector& aMom,const Material& aMaterial,const bool& aIsMultiScatt=false)
		{
			//need aMaterial.TR and aMaterial.L
			//aPos and aMom are inputs, also outputs
			//aPos: position aMom: momentum
			double lE=aMom(3);
			double ms_phi,ms_theta;
			if ( aIsMultiScatt && fabs(aMaterial.TR)>1e-15 )
			{
				ms_phi=MultiScattering(lE,aMaterial.TR);//rad
				ms_theta=MultiScattering(lE,aMaterial.TR);//rad
			}
			else {
				ms_phi=0;
				ms_theta=0;
			}
			//ms_theta=2*PI*gRandom->Rndm();
			double lth=aMom.X()/aMom.Z();
			double lph=aMom.Y()/aMom.Z();

			//pass L/2
			aPos.SetX(aPos.X()+aMaterial.L/2*lth);
			aPos.SetY(aPos.Y()+aMaterial.L/2*lph);
			aPos.SetZ(aPos.Z()+aMaterial.L/2);

			//change the angle and pass the rest L/2
			lth=(tan(ms_theta)+lth)/(1-tan(ms_theta)*lth);
			lph=(tan(ms_phi)+lph)/(1-tan(ms_phi)*lph);
			aMom.SetX(lth);
			aMom.SetY(lph);
			aMom.SetZ(1);
			aPos.SetX(aPos.X()+aMaterial.L/2*lth);
			aPos.SetY(aPos.Y()+aMaterial.L/2*lph);
			aPos.SetZ(aPos.Z()+aMaterial.L/2);
		}
		/*}}}*/

		/*double MultiScattering(const double& aE,const double& aTR){{{*/
		double MultiScattering(const double& aE,const double& aTR)
		{
			//only for electron
			double lPsq=aE*aE-ELECTRON_MASS*ELECTRON_MASS;
			double bcp=lPsq/aE;
			double ltheta0=13.6/bcp*sqrt(aTR)*(1+0.038*log(aTR));
			if ( aTR!=0 )
			{
				//return gRandom->Gaus(0,ltheta0/2.3548);//rad sigma=width/(2*sqrt(2)
				//ltheta/w, sg=sigma of y_tg in data
				//w=1.461e6*sg^2-6976*sg+9.316
				return gRandom->Gaus(0,ltheta0/1.3548);//rad sigma=width because of the
				//calculation of energy loss is behind the multiscattering
				//I add more multiscattering
			}
			else
				return 0;
		}
		/*}}}*/

		/*void SetMaterial(Material& aMaterial){{{*/
		void SetMaterial(Material& aMaterial)
		{
			aMaterial.M=aMaterial.A*AMU;//MeV
			if ( aMaterial.L==0 && aMaterial.rho!=0 ) {
				aMaterial.L=aMaterial.T/aMaterial.rho;
			}
			aMaterial.X0=Rad_Len(aMaterial.Z,aMaterial.A);
			if ( aMaterial.X0!=0 )
				aMaterial.TR=aMaterial.T/aMaterial.X0;
			else
				aMaterial.TR=0;
			aMaterial.bt=b(aMaterial.Z)*aMaterial.TR;
		}
		/*}}}*/

		/*Material GetMixture(const vector<Material>& aWin){{{*/
		Material GetMixture(const vector<Material>& aWin)
		{
			size_t i;
			size_t imax=aWin.size();
			Material mixture;
			mixture.Name="mixture";
			mixture.L=0;
			mixture.A=0;
			//get mixture TR
			for ( i=0; i<imax; i++ )
			{
				mixture.L+=aWin[i].L;
				mixture.A+=aWin[i].A;
			}
			mixture.TR=0;
			for ( i=0; i<imax; i++ )
			{
				mixture.TR+=aWin[i].A/mixture.A*aWin[i].TR;
			}
			return mixture;
		}
		/*}}}*/

		/*void GetRef_Plane(TLorentzVector& aoPos,TLorentzVector& aoMom,const vector<Material>& aWinBefore,const vector<Material>& aWinAfter,const double& aL,const double& aOffset){{{*/
		void GetRef_Plane(TLorentzVector& aoPos,TLorentzVector& aoMom,const vector<Material>& aWinBefore,const vector<Material>& aWinAfter,const double& aL,const double& aOffset)
		{
			//Get Refinement aoPos and aoMom for each plane on q1,q2,d,q3,fp
			//ao means input and output
			//aL=Vacuum Length, aOffset=distance between interaction point and z=0 in TCS
			vector<Material> AllWins;
			Material Vacuum;

			AllWins.clear();
			AllWins=aWinBefore;
			Vacuum.Name="Vacuum";
			Vacuum.L=aL;
			Vacuum.A=0;
			Vacuum.TR=0;
			AllWins.push_back(Vacuum);
			size_t i;
			for ( i=0; i<aWinAfter.size(); i++ )
				AllWins.push_back(aWinAfter[i]);
			Vacuum.L=0;
			//Printf("Pos(%g,%g,%g),(th=%g,ph=%g)",aoPos(0),aoPos(1),aoPos(2),aoMom(0)/aoMom(2),aoMom(1)/aoMom(2));
			for ( i = 0; i < AllWins.size(); ++i ) {
				Transport(aoPos,aoMom,AllWins[i],IsMultiScat);//move to mag
				//Printf("Pos(%g,%g,%g),(th=%g,ph=%g)",aoPos(0),aoPos(1),aoPos(2),aoMom(0)/aoMom(2),aoMom(1)/aoMom(2));
				Vacuum.L-=AllWins[i].L;
			}
			Vacuum.L+=aOffset;
			Transport(aoPos,aoMom,Vacuum,IsMultiScat);//back to ztg=0 in TCS
			//Printf("Pos(%g,%g,%g),(th=%g,ph=%g)",aoPos(0),aoPos(1),aoPos(2),aoMom(0)/aoMom(2),aoMom(1)/aoMom(2));
		}
		/*}}}*/

		/*double sigma_M(const double& aE,const double& aTheta){{{*/
		double sigma_M(const double& aE,const double& aTheta)
		{
			//Mott Cross Section
			//aE=MeV,aTheta=deg
			double ltheta=aTheta/2.*PI/180;
			double mott;
			if ( ltheta!=0 )
				mott=pow(ALPHA*cos(ltheta)/(2*aE*pow(sin(ltheta),2)),2);
			else
				mott=0;
			return mott*MEV2SR_TO_NBARNSR; //nbarn
			//return mott; //nbarn
		}
		/*}}}*/
  
		/*double CalcRValue(const double& ath,const double& aph,const double& ay,const double& adp){{{*/
		double CalcRValue(const double& ath,const double& aph,const double& ay,const double& adp)
		{
			int i,j,k;

			i=0;
			double* lxy=new double[NELEMENTS];
			lxy[i++]=ath;
			lxy[i++]=aph;
			lxy[i++]=ay;
			lxy[i++]=adp;
			double prod=-1000;
			if ( fNCuts>0 ) {
				double lnormal;
				double lomega;
				for ( i = 0; i < fNCuts; ++i ) {
					//slope*x+intersection-y=0
					
					lnormal=sqrt(fLineProperty[i][LINE_SLOPE]*fLineProperty[i][LINE_SLOPE]+1);//Get normalized factor (sqrt(a*a+b*b)) a=slope b=-1
					lomega=0;
					k=0;
					for ( j = 0; j < NELEMENTS; ++j ) {
						lomega+=fXY[i][j]*pow(-1.0,fXY[i][j]+k)*pow(fLineProperty[i][LINE_SLOPE],k)*lxy[j];
						//if ( i==1 ) {
						//	Printf("i=%d,j=%d,k=%d,lomega=%g",i,j,k,lomega);
						//	Printf("fXY[%d][%d]=%d,k=%d,fLineProperty[%d][%d]=%g,lxy[%d]=%g,lomega=%d*pow(-1,%d)*pow(%g,%d)*%g=%g",i,j,fXY[i][j],k,i,LINE_SLOPE,fLineProperty[i][LINE_SLOPE],j,lxy[j],fXY[i][j],fXY[i][j]+k,fLineProperty[i][LINE_SLOPE],k,lxy[j],fXY[i][j]*pow(-1.0,fXY[i][j]+k)*pow(fLineProperty[i][LINE_SLOPE],k)*lxy[j]);

						//}
						k+=fXY[i][j];
					}
					lomega+=fLineProperty[i][LINE_INTERSECTION];
					lomega/=lnormal*fLineProperty[i][LINE_SIGN];
					//Printf("lomega[%d]=%g",i,lomega);
					if ( i==0 ) {
						prod=lomega;
					}
					else {
						prod=PROD_AND(prod,lomega);
					}
				}
			}
			
			delete [] lxy;
			return prod;
		}
		/*}}}*/

		/*double PROD_AND(const double& ax,const double& ay){{{*/
		double PROD_AND(const double& ax,const double& ay)
		{
			double prod;
			prod=TMath::Min(ax,ay);
			//prod=ax+ay-sqrt(ax*ax+ay*ay);
			return prod;
		}
		/*}}}*/


	public:
		//Not safe, but it's ok for physicsists
		/*Member Data from file{{{*/
		int Id; //Event Id
		int IsPassed; //if pass through the magnet, 1=true, 0=false(not pass)
		int IsQualified; //if rec var is in the range of gen, 1=true, 0=false(not qualified)
		double theta; //scattering angle(deg)
		Material Target; //target material
		Material Win_i,Win_f;//front[inital] window and back[final] window, stick with target
		vector<Material> Win_Before_Mag;//windows after target and before magnetic
		vector<Material> Win_After_Mag;//windows after target and after magnetic
		double T_theta; //target angle(deg), the central line of target is along beam
		//----
		//|  |
		//|  |--->central line, so the T_theta is defined like phi_tg
		//|  |
		//----
		/*}}}*/

		/*Member Data derived from variables above{{{*/
		//TLorentzVector ( vec, t ) = ( k,E ) in ROOT
		//HCS: Hall Coordinate System
		//TCS: Target Coordinate System
		TLorentzVector s; //s(s3,E_s) 4-momentum position of the incident electron in HCS
		TLorentzVector target_edgepoint_TRCS; //4-momentum position of the edge point in HCS (if theta>0,T_H/2,if <0,-T_H/2)
		TLorentzVector s_TCS; //s_TCS(s3,E_s) 4-momentum position of the incident electron in TCS, and ztg=0
		TLorentzVector p_TCS; //p_TCS(p3,E_p) 4-momentum position of the outgoing electron in TCS, and ztg=0
		TLorentzVector p_inter_point_TCS; //p_inter_point_TCS(p3,E_p) 4-momentum position of the outgoing electron in TCS, and ztg=reactz_TCS
		TLorentzVector p_P; //p_P 4-momentum momentum of the outgoing electron in HCS
		TLorentzVector p_P_TCS; //p_P_TCS 4-momentum momentum of the outgoing electron in TCS

		double Angle; //real scattering angel, not from file,rad
		double Angle_Deg; //real scattering angel, not from file,Deg
		double sinsq_Angle;//sin(Angle/2)^2

		double theta_rad;//scattering angle in SAMC
		double sinsq;//sin(theta/2)^2
		double reactz_TCS;//
		double btr;//b*tr (equivalent radiator)

		double Q2;// MeV*MeV
		double q2;//q2=-Q2 MeV*MeV
		/*}}}*/

		/*Member Data from function{{{*/
		/*}}}*/

		/*Member Data from random{{{*/
		double beam_x; //cm
		double beam_y; //cm
		double reactz_gen; // interaction z position cm
		double E_s;//=incident beam energy (MeV)
		double E_p;//=incident beam energy (MeV)
		//gen=generator
		double th_tg_gen;//theta target
		double ph_tg_gen;//phi target
		double dp_gen;//dp at target
		/*}}}*/

		/*Member Data derived from variables above{{{*/
		double x_tg_gen; //
		double y_tg_gen; //
		/*}}}*/

		/*Member Data for refinement{{{*/
		double dp_ref; //dp for John.LeRose matrix
		double x_tg_ref; //x tg for John.LeRose matrix
		double y_tg_ref; //y tg for John.LeRose matrix
		double th_tg_ref; //th tg for John.LeRose matrix
		double ph_tg_ref; //ph tg for John.LeRose matrix
		TLorentzVector p_TCS_ref;//
		TLorentzVector p_P_TCS_ref;//
		/*}}}*/

		/*Member Data for focal plane {{{*/
		TLorentzVector p_FP; //p_FP 4-momentum momentum of the outgoing electron in focal plane
		TLorentzVector p_P_FP; //p_P_FP 4-momentum momentum of the outgoing electron in focal plane
		double x_fp; //
		double y_fp; //
		double th_fp; //
		double ph_fp; //
		double th_fp_no_ort; //th_fp without TXFIT orthogonalization
		double q1ex[2];//0:x 1:y
		double dent[2];
		double dext[2];
		double q3en[2];
		double q3ex[2];
		bool IsPassedQ1Ex;
		bool IsPassedDipoleEn;
		bool IsPassedDipoleEx;
		bool IsPassedQ3En;
		bool IsPassedQ3Ex;
		/*}}}*/

		/*Member Data for reconstruct target varibles {{{*/
		double x_tg_rec; //cm
		double y_tg_rec; //cm
		double th_tg_rec; //
		double ph_tg_rec; //
		double dp_rec;
		double reactz_rec;
		double rvalue;
		double cs_M;//mott cross section
		double cs_Final;//final cross section; 
		double Angle_rec;//Calculated from reconstructed variables,
		double Qsq; //Calculated from reconstructed variables, Q2 is from the generated variables 
		double Xbj;//Calculated from reconstructed variables,
		/*}}}*/

		/*void Init(){{{*/
		void Init()
		{
			IsPassed=1;
			IsQualified=1;
			Win_Before_Mag.clear();
			Win_After_Mag.clear();
		};
		/*}}}*/

		/*void Print(){{{*/
		void Print()
		{
			printf("--------------------\n");
			//refer to http://www.fileformat.info/info/unicode/char/search.htm for \x
			printf("%*s=(%*.2f,%*.2f,%*.4f,%*.2f) (4-momentum position of target edgepoint [HCS])\n",20,"target_edgepoint_TRCS(x,y,z,E_s)",10,target_edgepoint_TRCS(0),10,target_edgepoint_TRCS(1),10,target_edgepoint_TRCS(2),10,target_edgepoint_TRCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of incident electron [HCS])\n",20,"s(x,y,z,E_s)",10,s(0),10,s(1),10,s(2),10,s(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of incident electron [TCS],ztg=0)\n",20,"s_TCS(x,y,z,E_s)",10,s_TCS(0),10,s_TCS(1),10,s_TCS(2),10,s_TCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of outgoing electron [TCS], ztg=0)\n",20,"p_TCS(x,y,z,E_s)",10,p_TCS(0),10,p_TCS(1),10,p_TCS(2),10,p_TCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of outgoing electron [TCS], ztg=reactz_TCS=%.2f)\n",20,"p_inter_point_TCS(x,y,z,E_s)",10,p_inter_point_TCS(0),10,p_inter_point_TCS(1),10,p_inter_point_TCS(2),10,p_inter_point_TCS(3),reactz_TCS);
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [HCS])\n",20,"p_P(x,y,z,E_s)",10,p_P(0),10,p_P(1),10,p_P(2),10,p_P(3));
			printf("%*s=(%*.2e,%*.2e,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [TCS])\n",20,"p_P_TCS(x,y,z,E_s)",10,p_P_TCS(0),10,p_P_TCS(1),10,p_P_TCS(2),10,p_P_TCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of outgoing electron [TCS], ztg=0 after refinement)\n",20,"p_TCS_ref(x,y,z,E_s)",10,p_TCS_ref(0),10,p_TCS_ref(1),10,p_TCS_ref(2),10,p_TCS_ref(3));
			printf("%*s=(%*.2e,%*.2e,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [TCS] after refinement)\n",20,"p_P_TCS_ref(x,y,z,E_s)",10,p_P_TCS_ref(0),10,p_P_TCS_ref(1),10,p_P_TCS_ref(2),10,p_P_TCS_ref(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [FP])\n",20,"p_P(x,y,z,E_s)",10,p_P(0),10,p_P(1),10,p_P(2),10,p_P(3));
			printf("%*s=(%*.2e,%*.2e,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [FP])\n",20,"p_P_FP(x,y,z,E_s)",10,p_P_FP(0),10,p_P_FP(1),10,p_P_FP(2),10,p_P_FP(3));
			printf("%-*s=%*d %-*s %-*s\n",  15,"Id",        10,Id,        8,"",              40,"(Event Id)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"E_s",       10,E_s,       8,"MeV",           40,"(Es Incident Energy)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"E_p",       10,E_p,       8,"MeV",           40,"(Ep Scattering Energy)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"theta",     10,theta,     9,"\xc2\xb0",      40,"(\xce\xb8 scattering angle deg)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"theta_rad", 10,theta_rad, 8,"",              40,"(\xce\xb8 scattering angle rad)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"T_theta",   10,T_theta,   9,"\xc2\xb0",      40,"(\xce\xb8t target angle)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"sinsq",     10,sinsq,     8,"",              40,"(sin\xc2\xb2(\xce\xb8/2))");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"M_T",       10,Target.M,  8,"MeV",           40,"(Mass of target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"beam_x",    10,beam_x,    8,"cm",            40,"(Beam X)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"beam_y",    10,beam_y,    8,"cm",            40,"(Beam Y)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"reactz_gen",    10,reactz_gen,    8,"cm",            40,"(Z react)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"x_tg_gen",      10,x_tg_gen,      8,"cm",            40,"(x target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"y_tg_gen",      10,y_tg_gen,      8,"cm",            40,"(y target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"th_tg_gen",     10,th_tg_gen,     8,"",              40,"(tan(\xce\xb8tg))");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"ph_tg_gen",     10,ph_tg_gen,     8,"",              40,"(tan(\xcf\x86tg))");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"dp_gen",        10,dp_gen*100,    8,"%",             40,"(dp at target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"th_tg_ref", 10,th_tg_ref, 8,"",              40,"(tan(\xce\xb8tg))[ref_tg]");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"ph_tg_ref", 10,ph_tg_ref, 8,"",              40,"(tan(\xcf\x86tg))[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"dp_ref",    10,dp_ref*100,8,"%",             40,"(dp at target)[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"x_tg_ref",  10,x_tg_ref,  8,"cm",            40,"(x target)[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"y_tg_ref",  10,y_tg_ref,  8,"cm",            40,"(y target)[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"x_fp",  10,x_fp,  8,"cm",            40,"(x target)[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"y_fp",  10,y_fp,  8,"cm",            40,"(y target)[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"th_fp", 10,th_fp, 8,"",              40,"(tan(\xce\xb8tg))[ref_tg]");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"ph_fp", 10,ph_fp, 8,"",              40,"(tan(\xcf\x86tg))[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"x_tg_rec",      10,x_tg_rec,      8,"cm",            40,"(x target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"y_tg_rec",      10,y_tg_rec,      8,"cm",            40,"(y target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"th_tg_rec",     10,th_tg_rec,     8,"",              40,"(tan(\xce\xb8tg))");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"ph_tg_rec",     10,ph_tg_rec,     8,"",              40,"(tan(\xcf\x86tg))");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"dp_rec",        10,dp_rec*100,    8,"%",             40,"(dp at target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"rvalue",        10,rvalue,    8,"",             40,"(rvalue of target reconstructed variables)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Angle",     10,Angle,     8,"rad",           40,"(real scattering angle[rad])");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Angle",     10,Angle*RadToDeg(),     9,"\xc2\xb0",           40,"(real scattering angle[deg])");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"Q2",        10,Q2,        9,"MeV\xc2\xb2",   40,"(Q\xc2\xb2)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"q2",        10,q2,        9,"MeV\xc2\xb2",   40,"(-Q\xc2\xb2)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"btr",       10,btr,       8,"rad_len",       40,"(b*tr equivalent radiator unit in rad_len)");
			PrintMaterial(Target);
			PrintMaterial(Win_i);
			PrintMaterial(Win_f);
			size_t i;
			printf("\t\tWindows List Before Magnetic\n");
			for ( i=0; i<Win_Before_Mag.size(); i++ )
			{
				PrintMaterial(Win_Before_Mag[i]);
			}
			printf("\t\tWindows List After Magnetic\n");
			for ( i=0; i<Win_After_Mag.size(); i++ )
			{
				PrintMaterial(Win_After_Mag[i]);
			}

		};
		/*}}}*/
		
		/*void PrintMaterial(){{{*/
		void PrintMaterial(const Material& aMaterial)
		{
			printf("%-*s=%*d %-*s %-*s\n",15,Form("%s_Z",aMaterial.Name.c_str()),     10,aMaterial.Z,       8,"",  40,"(Atomic Number)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,Form("%s_A",aMaterial.Name.c_str()),     10,aMaterial.A,       8,"g/mol",  40,"(Atomic Weight)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,Form("%s_L",aMaterial.Name.c_str()),     10,aMaterial.L,       8,"cm",  40,"(length)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,Form("%s_T",aMaterial.Name.c_str()),     10,aMaterial.T,       9,"g/cm\xc2\xb2",  40,"(thickness)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,Form("%s_TR",aMaterial.Name.c_str()),     10,aMaterial.TR,       8,"rad_len",  40,"(thickness in rad_len)");
			printf("%-*s  =%*.2e %-*s %-*s\n",16,Form("%s_\xf0\x9d\x9c\x8c",aMaterial.Name.c_str()),     10,aMaterial.rho,       9,"g/cm\xc2\xb3",  40,"(density)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,Form("%s_X0",aMaterial.Name.c_str()),     10,aMaterial.X0,       9,"g/cm\xc2\xb2",  40,"(Radiation Length)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,Form("%s_bt",aMaterial.Name.c_str()),     10,aMaterial.bt,       8,"",  40,"(bt)");
		}
		;/*}}}*/
		
}
;/*}}}*/
