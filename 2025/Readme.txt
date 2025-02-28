%%20250228 업데이트 - 김준연

원본 데이터 위치 : 'G:\공유 드라이브\BSL-Data\Data\Hyundai_dataset\현대차파우치셀 (rOCV,Crate)\NE_MCT25oC_HPPC25oC_OCV_KENTECH_송부.xlsx' (현대차 제공)

코드 목적 (연구 목적) 

	기존 RPT(원리실험, HPPC)에서 전압 거동을 fitting을 통해 시상수를 추출하는 방식은 많이 진행되었음
	여기서 추출한 RC을 validation으로 주행 데이터에 적용하여 검증하는 방향성은 많이 제시되었지만, 즉 기존 평가 실험은 실시간 배터리 상태 추정에는 적합하지 않음 
	--> 실시간으로 주행 데이터로부터 얻은 시상수가 어느정도 합리적인 시상수라면, 실시간 진단이 가능할것이라 예상됨 
	--> 즉, 실시간 주행으로 부터 얻은 time scale이 HPPC에서 얻은 time scale과 비교하여 완벽하게 맞지는 않지만 어느정도 겹치는 시상수가 존재하는것을 증명하는것이 목적
	--> 물론 비슷한 SOC에서 얻은 저항거동으로 비교분석해야함 (e.g HPPC 저항 (SOC 92) ~~ Driving 저항 (SOC 93) 비교)

데이터 분석

	기존 현대차 제공 데이터 : MCT-1,2,3,4,5,6 + HPPC + OCV

	MCT data  : 주행 data
		1,2,4,5 : 도로 주행 데이터 ( 도심 주행 + 고속 섞임)
		3,6 : 고속 주행 데이터 (108km/h로 등속 주행) --> CC 방전이라 볼수 있으므로 RC fitting은 진행하지 않음
		1RC, 2RC 코드 완료 ( 추후 DRT(nRC) , RC - W (와버그) 코드 추가 필요)
		Vmodel = OCV(SOC) + IR0 + Vrc1 + Vrc2 ( 주의 : discrete한 식, dt 사용)
		정확한 모델 추정을 위해서는 OCV에 대한 오차를 줄이는것이 중요 
		--> BMS SOC(원본 엑셀 데이터) or  재계산 SOC(Rest 전압으로부터 역추적한 SOC0 + 쿨롱 카운팅) 선택 필요
		Fitting output : [R0,R1,R2,Tau1, Tau2] (2RC 기준)


	HPPC data : RPT data
		저항 거동 분석 
		SOC 이동용 전류와 저항거동 분석용 pulse을 구별하는것이 중요 (저항 분석은 pulse에서만 분석) 
		---> 어떻게 pulse을 주고, 어떻게 soc을 이동했는지 분석한 page : 														(https://docs.google.com/presentation/d/1cnDEUu_Vttk9ZGytiZJpc0Oq5dJLfzs16MD2vxIUBAc/edit#slide=id.g32ff5fc5509_3_0)
		Fitting에 사용될 Pulse 구조체 (MCT soc에 대응되는 pulse) 을 생성 , 펄스 시작, 끝 Vrest을 찾아서 대응되는 SOC 양끝을 point out 한후, 양 끝점을 이음 ( 주의 : 분모 Q_batt가 아닌 Q_trip 사용해야함)
		Delta V = IR0 + IR1(1-exp(t/tau1)) + IR2(1-exp(t/tau2)) (주의 : continous 한 식, 누적식임)
		Fitting output : [R0, R1, R2, Tau1, Tau2] (2RC 기준)
		
	OCV data : RPT data
		SOC-OCV 곡선 얻을 수 있음
		cell 곡선, pack 곡선 존재
		Pack은 192s 2p로 이루어져 있음


---> 2가지의 fitting 결과, 특히 Tau1, Tau2을 비교 분석하는것이 목적
---> fitting 진행시, 초기값에 영향을 많이 받는것이 밝혀짐 , 아마 추후 toleracnce 을 조정하여 fitting robust 필요 예상 
---> 재계산 SOC (쿨롱 카운팅 기반) 을 사용했을때, 시작점의 정확성 때문에 BMS SOC보다 더 정확한 fitting 결과값을 나타냄 


자세한 코드 설명

초기값 설정 ref

	R0 : 1m옴 
	R1 : 0.5m옴
	R2 : 0.5m옴 (기존 교수님이 HPPC RC fitting 결과값 기준 사용)
	Tau1 : 짧은 scale (전압 거동이 초반에 팍 튐, 짤은 scale을 포착하기 위함), (CC_1RC,2RC_cost에서 초기값 가져옴)
	Tau2 : 긴 scale (뒤쪽 전압 거동이 CT에 의해 앞쪽 전압거동과 다른 거동을 보임, 상대적으로 더 긴 scale을 사용)  (CC_1RC,2RC_cost에서 초기값 가져옴)

0. rOCV hyundai.m

input : 현대차에서 제공한 RPT data 
output : Q_batt을 rOCV에서 CC계산하여 가져옴 (Q_batt : 56.2396) %[Ah]

1. ECM-MCT (BMS SOC만 그림), ECM-MCT2 (BMS SOC + 재계산 SOC 둘다 그림)

input : 엑셀 파일 (현대차 제공)
output : 가공된 MCT data (MCT_Results)

데이터 분석 목적

2. MCT_fitting_CC_1RC,MCT_fitting_CC_2RC

input : 가공된 MCT data
output : fitting 결과값, fitting plot

재계산 soc 기준으로 fitting한 code

3. MCT_fitting_BMS_1RC,MCT_fitting_CC_2RC

input : 가공된 MCT data
output : fitting 결과값, fitting plot

원본 BMS SOC 기준으로 fitting 한 code 

4. CC_1RC_cost, CC-2RC_cost, BMS_1RC_cost

input : 가공된 MCT data
output : RMSE vs Tau1, Tau2 

fitting 결과값이 초기값의 영향을 많이 받기 때문에 fitting이 정말 잘 되었는지 확인할 필요가 있음
(현재 code에서는 초기값을 지정하면 fitting 후 결과값이 그 값을 계속 따라가는 경향성이 있음)

원래 fitting 해야할 parameter는 [R0,R1,R2,Tau1,Tau2] 지만, Tau1과 Tau2을 grid로 범위를 주어서 (tau1의 범위 : 0.1-20, Tau2의 범위 : 10 - 180초)
tau1, tau2 조합으로 R0, R1, R2 fitting 하고, 반복 시행 후  COST(tau1, tau2) surface 그리기 --> 최적의 Tau1, Tau2을 초기값으로 사용하기 










	
	
	