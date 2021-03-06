

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
%ctrPacketLength = 200;
diffPacketLength = 50; % 차분 패킷
ctrPacketLength = 200;
%diffPacketLength = 1; % 차분 패킷
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
rmax=900;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%Computation of do

do=sqrt(Efs/Emp);

% 초기값 처리 여부
IS_INITIL_LEACH = true;

%병합 처리 여부
IS_MERGE = true;

cluster_data_count = 20;
leach_data = [];
initil_leach_data = [];
leach_data_length = zeros(1, rmax);
leach_node_data_length = zeros(1, rmax);
initil_data_length = zeros(1, rmax);
initil_node_data_length = zeros(1, rmax);
pro_data_length = zeros(1, rmax);
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

for leach_round=1:1:2
    for i=1:1:n
        %initially there are no cluster heads only nodes
        S(i).type='N';
        S(i).E=Eo;
        S(i).ENERGY=0;
        S(i).G=0;
        S(i).cl=0;
        dead_node_id = zeros(1, n);
    end
    S(n+1).xd=sink.x;
    S(n+1).yd=sink.y;
    if leach_round ~= 1
        IS_INITIL_LEACH = false;
       
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
    flag_hafe_dead=0; 
    flag_all_dead=0; 
    
    
    r=-1;
    while r<=rmax
        r = r+1;
     % r
      % make packetData for 20171101 minute data
      round_sensing_data = [];
      unzip_round_sensing_data = [];
      round_data = zeros(1, cluster_data_count);
      live_node = 0;
      
      %for i=1:1:n
      % if (S(i).E>0) 
      %  live_node = live_node +1;
      % end
      %end
      %c_node = 1;
      %if live_node * p > 0
      %    c_node = live_node * p;
      %end
      %cluster_data_count = int8(live_node / c_node);
      
      if r ~= 0
          diff_row_val = (sensing_data(r+1, 1:cluster_data_count)*10)-(sensing_data(r, 1:cluster_data_count)*10);
      end 
      for i=1:1:cluster_data_count
         if (IS_INITIL_LEACH )
            if r == 0
                % 차분값
                round_sensing_data = [round_sensing_data ,(i), (sensing_data(r+1, i)*10)];
            else
              if (IS_MERGE)
                  % 병합처리시 변화 없는 데이터 전송 안함
                  if diff_row_val(i) ~= 0
                    round_sensing_data = [round_sensing_data ,(i), (diff_row_val(i))];
                  % round_data(i) = diff_row_val(i);
                  end
              else
                  round_sensing_data = [round_sensing_data ,(i), (diff_row_val(i))];
              end   
            end
         else 
             % 기본 leach의 경우
            round_sensing_data = [round_sensing_data (i) (sensing_data(r+1, i)*10)];
         end
        % (sensing_data(1, :)*10 ) - (sensing_data(2, :)*10 )
        % sensing_data([1,2], r:i)*10
      end
      unzip_round_sensing_data = round_sensing_data;
      

      %Operation for epoch
      can_be_cluster_header_cnt = 0;
      live_node_cnt = 0;
      for i=1:1:n
          if S(i).G <= 0
              can_be_cluster_header_cnt = can_be_cluster_header_cnt + 1;
          end
          if S(i).E>0
              live_node_cnt = live_node_cnt +1;
          end
      end
      if can_be_cluster_header_cnt == 0
          for i=1:1:n
            S(i).G=0;
            S(i).cl=0;
         end
      end
      
      
      
      if r > 880 && r < 890
          for i=1:1:n
              if S(i).E > 0
                   S(i).E;
              end
          end
      end
      
      if(mod(r, round(1/p))==0)
      %if(mod(r, 21)==0)
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
    % 병합 처리시 죽은 노드의 특수 ID 값 추가 
    if ( r ~= 0 && IS_INITIL_LEACH && IS_MERGE )
    %if false
       round_dead_node_id = [];
       this_round_dead_node_id = dead_node_id == 1;
       for i=1:1:n
           if this_round_dead_node_id(i) == 1
               round_dead_node_id = [round_dead_node_id i];
           end
       end
       % 죽은 노드 ID 저장
       if ~strcmp(round_dead_node_id, '')
        round_sensing_data = [round_sensing_data round_dead_node_id];
       end
       % 데이터 병합값 생성
       % 기준 값
       %round_data_min = min(round_data);
       % 기준값으로 차분값 처리
       %round_data = round_data-round_data_min;
       %round_data_unique = unique(round_data);
       
       % 병합처리
       % strrep(strtrim(cellstr(num2str(A))), '  ', '|')
       
       %for i=1:1:length(round_data_unique)
       %    ids_val =  char(regexprep(strtrim(cellstr(num2str(find(round_data==round_data_unique(i))))), '\s+', '|'));
           
       %    if round_data_unique(i) == 0
       %        round_sensing_data = sprintf('%s0|%s:%d,',round_sensing_data, ids_val, round_data_min);
       %    else
              
       %        if (round_data_unique(i) - abs(round_data_min)) ~= 0
       %         round_sensing_data = sprintf('%s%s:%d,',round_sensing_data, ids_val, round_data_unique(i));
       %        end

       %   end
       % end
       round_sensing_data;
       % round_sensing_data = sprintf('%s%d:%d,',round_sensing_data ,(i*10), (diff_row_val(i)));
       
    end
    
    
    % length(dec2bin(round_sensing_data, 16) - '0')
    % string to 16 bits
    
    % packetLength = length(dec2bin(round_sensing_data, 16) - '0')*16;
    [round_row,round_column] = size(dec2bin(abs(round_sensing_data)));
    round_column = round_column+1; % 음수 처리
    packetLength = round_row * round_column;
    
    
    %packetLength = packetLength * 2;
    %if leach_round == 1
    %    packetLength = 100000
    %else 
    %    packetLength = 50000
    %end
   %if IS_MERGE  && leach_round == 2
    % packetLength = compressionLZW(round_sensing_data);
   %end
    % 0~20 까지 데이터 출력용
    if r < 20
        %fprintf('round :%d\n', leach_round); 
        %fprintf('packetLen:%d, data: ',packetLength); 
        %round_sensing_data
    end
    
    
    if leach_round == 1 
        leach_data_length((r+1)) = packetLength;
    elseif leach_round == 2
        pro_data_length((r+1)) = packetLength;
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
            fprintf('dead first :%d\n', r);
        end
 elseif (dead == 100)
         if(flag_hafe_dead==0)
             flag_hafe_dead=1;
            fprintf('dead hafe :%d\n', r); 
         end
    elseif (dead == n)
         if(flag_all_dead==0)
             flag_all_dead=1;
            fprintf('dead all :%d\n', r); 
         end        
    end

    countCHs=0;
    cluster=1;
    for i=1:1:n
       if(S(i).E>0)
         temp_rand=rand;     
         if ((S(i).G)<=0)
         %if true
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
                if ( r ~= 0 && IS_INITIL_LEACH && IS_MERGE && leach_round == 3)
                    S(i).E = S(i).E-(EDA * (packetLength));
                end
                
                if r == 0
                     S(i).E 
                 end
                packets_TO_BS = packets_TO_BS+1;
                PACKETS_TO_BS(r+1) = packets_TO_BS;
            end     
         end
       end 
    end

    if (live_node_cnt > 0) && (countCHs == 0)
        r = r-1;
        for i=1:1:n
            S(i).G=0;
            S(i).cl=0;
         end
        continue;
    end
    
    STATISTICS(r+1).CLUSTERHEADS = cluster-1;%
    CLUSTERHS(r+1)= cluster-1;


    
    
    %Election of Associated Cluster Head for Normal Nodes
    for i=1:1:n
       if (S(i).type=='N' && S(i).E>0) 
        % min_dis = sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%
         min_dis = INFINITY; 
         if (cluster-1>=1)
         %if true
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
                if r == 0
                    nodeData =  [i (sensing_data(r+1, i)*10)];
                else
                    nodeData = [i ((sensing_data(r+1, i)*10)-(sensing_data(r, i)*10))];
                end
            else
                nodeData =[i (sensing_data(r+1, i)*10)];
            end

            nodePacketLength = numel(dec2bin(abs(nodeData)));
            
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
             if ( r ~= 0 && IS_INITIL_LEACH && IS_MERGE )
                 S(i).E = S(i).E-(EDA*(nodePacketLength*2)); % 차분 처리 에너지 소비 차분/갱신
             end
             if r == 0
                 S(i).E 
             end
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
         lzw_data = [x;y];
    elseif  leach_round == 2
        leach_data = [x;y];
    else
        initil_leach_data = [x;y];
    end
end

%plot(x,y,'r',x,z,'b');

l = plot(leach_data(1, [1:rmax]), leach_data(2, [1:rmax]), 'k--', lzw_data(1, [1:rmax]), lzw_data(2, [1:rmax]),'k-');
l(1).LineWidth = 1;
l(2).LineWidth = 1;

xlabel('Round');
ylabel('Number of Live Node');
legend('LEACH', 'Node ID Base DDP','Location','southwest');
hold on;

fix(mean(leach_data_length))
fix(mean(pro_data_length))
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