

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
sensing_data = csvread('day_data.csv');
%Field Dimensions - x and y maximum (in meters)
xm = 300;
ym = 300;
%x and y Coordinates of the Sink
sink.x =0.5 * xm;
sink.y = ym + 50;
%sink.x=50;
%sink.y=175;
%sink.x=0.5*xm;
%sink.y=0.5*ym;

%Number of Nodes in the field
n = 200;
%Optimal Election Probability of a node to become cluster head
p=0.05;
packetLength =6400;
ctrPacketLength = 200;
diffPacketLength = 50; % 차분 패킷
%Energy Model (all values in Joules)
%Initial Energy 
Eo = 0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

INFINITY = 999999999999999;
%maximum number of rounds
rmax=1000;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%Computation of do

do=sqrt(Efs/Emp);

% 초기값 처리 여부
IS_INITIL_LEACH = false;

%병합 처리 여부
IS_MERGE = true;

cluster_data_count = 20;
leach_data = [];
initil_leach_data = [];
leach_data_length = zeros(1, rmax);
initil_data_length = zeros(1, rmax);
lzw_data_length = zeros(1, rmax);
dead_node_id = zeros(1, n);


%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    % hold on;
end

for leach_round=1:1:3
    for i=1:1:n
        %initially there are no cluster heads only nodes
        S(i).type='N';
        S(i).E=Eo;
        S(i).ENERGY=0;
    end
    S(n+1).xd=sink.x;
    S(n+1).yd=sink.y;
    if leach_round ~= 1
        IS_INITIL_LEACH = true;
        dead_node_id = zeros(1, n);
    end
    %First Iteration
    figure(1);

    %counter for CHs
    countCHs=0;
    %counter for CHs per round
    rcountCHs=0;
    cluster=1;

    countCHs;
    rcountCHs=rcountCHs+countCHs;
    flag_first_dead=0; 
    
    
    
    for r=0:1:rmax 
     % r
      % make packetData for 20171101 minute data
      round_sensing_data=[];
      unzip_round_sensing_data=[];
      round_data = zeros(1, cluster_data_count);
      
      diff_cluster_val = sensing_data(r+1, cluster_data_count+4)*10;
      round_sensing_data = [diff_cluster_val];
      for i=1:1:cluster_data_count
          % cluster header sesing value
          
         if (IS_INITIL_LEACH )
            % 차분 처리
            diff_node_val = (sensing_data(r+1, i)*10)-diff_cluster_val;
            %if diff_node_val < 0
                round_sensing_data = [round_sensing_data diff_node_val];
            %else
            %    round_sensing_data = [round_sensing_data ((sensing_data(r+1, i)*10)-diff_cluster_val)];
            %end
         else 
             % 기본 leach의 경우
            round_sensing_data = [round_sensing_data (sensing_data(r+1, i)*10)];
         end
       end
      unzip_round_sensing_data = round_sensing_data;
      

      %Operation for epoch
      %can_be_cluster_header_cnt = 0;
      %for i=1:1:n
      %    if S(i).G == 0
      %        can_be_cluster_header_cnt = can_be_cluster_header_cnt + 1;
      %    end
      %end
      %if can_be_cluster_header_cnt == 0
      %    S(i).G=0;
      %    S(i).cl=0;
      %end
      %if(mod(r, round(1/p))==0)
      if(mod(r, 22)==0)
         for i=1:1:n
            S(i).G=0;
            S(i).cl=0;
         end
      end

    hold off;

    %Number of dead nodes
    dead=0;

    %counter for bit transmitted to Bases Station and to Cluster Heads
    packets_TO_BS=0;
    packets_TO_CH=0;
    %counter for bit transmitted to Bases Station and to Cluster Heads per round
    PACKETS_TO_CH(r+1)=0;
    PACKETS_TO_BS(r+1)=0;

    figure(1);
    
    
    for i=1:1:n
        %checking if there is a dead node
         if (S(i).E<=0)
           dead=dead+1;
           dead_node_id(i) =  dead_node_id(i)+1;
         end

         if (S(i).E>0)
            S(i).type='N';
         end
    end
     
    % length(dec2bin(round_sensing_data, 16) - '0')
    % string to 16 bits
    
    % packetLength = length(dec2bin(round_sensing_data, 16) - '0')*16;
    [round_row,round_column] = size(dec2bin(abs(round_sensing_data)));
    round_column = round_column+1; % 음수 처리
    packetLength = round_row * round_column;
    
    if leach_round == 3
     packetLength = numel(dec2bin(abs(round_sensing_data)));
    end
    if IS_MERGE  && leach_round == 2
        packetLength = compressionLZW(round_sensing_data);
    end
    % 0~20 까지 데이터 출력용
    if r < 20
        fprintf('round :%d\n', leach_round); 
        fprintf('packetLen:%d, data: ',packetLength); 
        jsonencode(round_sensing_data)
    end
    
    if leach_round == 1 
        leach_data_length((r+1)) = packetLength;
    elseif leach_round == 2
        lzw_data_length((r+1)) = packetLength;
    else
        initil_data_length((r+1)) = packetLength;
    end 
    
    %if (dead == n)
    %    r
    %end
    if r == rmax
       break;
    end

    STATISTICS(r+1).DEAD=dead;
    DEAD(r+1)=dead;

    %When the first node dies
    if (dead==1)
        if(flag_first_dead==0)
            first_dead=r;
            flag_first_dead=1;
        end
    end

    countCHs=0;
    cluster=1;
    for i=1:1:n
       if(S(i).E>0)
         temp_rand=rand;     
         if ((S(i).G)<=0) 
            %Election of Cluster Heads
            if(temp_rand <=(p/(1-p*mod(r,round(1/p)))))
                countCHs = countCHs+1;

                S(i).type = 'C';
                S(i).G = round(1/p)-1;
                C(cluster).xd = S(i).xd;
                C(cluster).yd = S(i).yd;

                distance=sqrt((S(i).xd-(S(n+1).xd))^2 + (S(i).yd-(S(n+1).yd))^2);%??sink?????

                C(cluster).distance = distance;
                C(cluster).id = i;
                X(cluster)=S(i).xd;
                Y(cluster)=S(i).yd;
                cluster=cluster+1;
                % 클러스터 헤드가 방송이 되었습니다.
                distanceBroad = sqrt(xm*xm+ym*ym);
                if (distanceBroad >=do)
                    S(i).E = S(i).E-(ETX*ctrPacketLength + Emp*ctrPacketLength*(distanceBroad*distanceBroad*distanceBroad*distanceBroad));%?????????
                else
                    S(i).E = S(i).E-(ETX*ctrPacketLength + Efs*ctrPacketLength*(distanceBroad*distanceBroad)); 
                end
                %Calculation of Energy dissipated 클러스터 헤드는 자체 패킷 에너지 소비량을 보냅니다.
                distance;
                % 병합 에너지 처리 제외
                if(distance>=do)
                     % S(i).E = S(i).E-((ETX+EDA)*packetLength+ Emp*packetLength*(distance*distance*distance*distance ));
                     S(i).E = S(i).E-((ETX)*packetLength+ Emp*packetLength*(distance*distance*distance*distance ));
                else
                     % S(i).E = S(i).E-((ETX+EDA)*packetLength+ Efs*packetLength*(distance*distance)); 
                     S(i).E = S(i).E-((ETX)*packetLength+ Efs*packetLength*(distance*distance)); 
                end
                % 병합 에너지 처리
                unzipPacketLength = numel(dec2bin(abs(unzip_round_sensing_data)));
                S(i).E = S(i).E-(EDA * unzipPacketLength);
                % S(i).E = S(i).E-(EDA * (packetLength)); % 차분 처리시
                % 차분 데이터 압축의 경우 패키지 사이즈로 에너시 소비 추가
                if ( r == 2  && IS_INITIL_LEACH && IS_MERGE )
                    S(i).E = S(i).E-(EDA * (packetLength));
                end
                
                
                packets_TO_BS = packets_TO_BS+1;
                PACKETS_TO_BS(r+1) = packets_TO_BS;
            end     
         end
       end 
    end

    STATISTICS(r+1).CLUSTERHEADS = cluster-1;%
    CLUSTERHS(r+1)= cluster-1;

    %Election of Associated Cluster Head for Normal Nodes
    for i=1:1:n
       if (S(i).type=='N' && S(i).E>0) 
        % min_dis = sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%
         min_dis = INFINITY; 
         if(cluster-1>=1)
             min_dis_cluster = 1;

             for c = 1:1:cluster-1 %
                %temp = min(min_dis,sqrt( (S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2 ) );
                temp = sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2);
                if (temp<min_dis)
                    min_dis = temp;
                    min_dis_cluster = c;
                end

                S(i).E = S(i).E - ETX * ctrPacketLength;
             end

             
             % Node data
            if (IS_INITIL_LEACH )
                nodeData = abs((sensing_data(r+1, i)*10)-sensing_data(r+1, cluster_data_count+4)*10);
            else
                nodeData = (sensing_data(r+1, i)*10);
            end

         
            nodePacketLength = numel(nodeData);
            
             %Energy dissipated by associated Cluster Head
             min_dis;
             if (min_dis > do)
                 S(i).E = S(i).E - (ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis * min_dis * min_dis * min_dis)); %
                 S(i).E = S(i).E - (ETX*(nodePacketLength) + Emp*nodePacketLength*( min_dis * min_dis * min_dis * min_dis)); %
             else
                S(i).E = S(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength*( min_dis * min_dis)); %
                S(i).E = S(i).E -(ETX*(nodePacketLength) + Efs*nodePacketLength*( min_dis * min_dis)); %
             end
             S(i).E = S(i).E - ETX*(ctrPacketLength);  %
             %if ( IS_INITIL_LEACH && IS_MERGE )
             %    S(i).E = S(i).E-(EDA*diffPacketLength); % 차분 처리 에너지 소비
             %end

             %Energy dissipated 
             if(min_dis > 0)
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA)*nodePacketLength ); % sensing energy
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ERX *ctrPacketLength ; %
                if (min_dis > do)%?
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis * min_dis * min_dis * min_dis));
                else
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( ETX*(ctrPacketLength) + Efs * ctrPacketLength*( min_dis * min_dis));
                end
               PACKETS_TO_CH(r+1) = n - dead - cluster + 1; 
             end

             S(i).min_dis = min_dis;
             S(i).min_dis_cluster = min_dis_cluster;

         end
      end
    end
    %hold on;

    countCHs;
    rcountCHs = rcountCHs + countCHs;
    figure(11)
    warning('OFF');
    [vx,vy]=voronoi(X(:),Y(:));
    plot(X,Y,'r+',vx,vy,'m-');
    hold on;
    voronoi(X,Y);
    axis([10 xm 0 ym]);


    end

    x=1:1:r;
    y=1:1:r;
    %z=1:1:r;

    for i=1:1:r
        x(i)=i;
        y(i) = n - STATISTICS(i).DEAD;
        %z(i)=CLUSTERHS(i);
    end
    if leach_round == 1
        leach_data = [x;y];
    elseif  leach_round == 2
        lzw_data = [x;y];
    else
        initil_leach_data = [x;y];
    end
end

%plot(x,y,'r',x,z,'b');

% plot(leach_data(1, [1:rmax]), leach_data(2, [1:rmax]), 'b:', lzw_data(1, [1:rmax]), lzw_data(2, [1:rmax]), 'g--', initil_leach_data(1, [1:rmax]), initil_leach_data(2, [1:rmax]), 'r-');
plot(leach_data(1, [1:rmax]), leach_data(2, [1:rmax]), 'b--', initil_leach_data(1, [1:rmax]), initil_leach_data(2, [1:rmax]), 'r-');
xlabel('Round');
ylabel('Number of Live Node');
legend('LEACH', 'CH Base DDP','Location','southwest');
%legend('Normal', 'LZW+Pro','Proposal');
hold on;

fix(mean(leach_data_length))
fix(mean(lzw_data_length))
fix(mean(initil_data_length))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS GRAPH PLOT SIR   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round
%  first_dead: the round where the first node died                   
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sink(50,175) ,ctrPacketLength=200,packetLength=4000,Eo=2J.