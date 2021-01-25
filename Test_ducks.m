clear all
%total duration of simulation
TMAX = 600;
%Change this number to change number of simulations to run
n_sims = 5;
plots = zeros(n_sims,4,TMAX); %array setup to store info about ducks each time step
%Starting simulation
for sims = 1:n_sims
fprintf('Simulation %d of %d \n', sims,n_sims)
%set agent speed/number 
v=0.3;
N = 100;
%Sets infection radius
r_transmission = 3.0;
duck_properties = zeros(N,13); %Duck property array: N x X length where N is number of ducks
                              %Duck properties being stored:
                              %1) x_pos, 2) y_pos, 3) random number for duck
                              %stop, 4) random number for duck landing, 5) angle
                              %of direction, 6) infected tracker, 7) time spent
                              %infected (0=no/1=yes), 8) immunity? (0=no/1=yes), 9) time spent
                              %paused, 10) time spent grounded,
                              %11) has the duck been infected this time
                              %step (0=no/1=yes) 12) has the duck been infected this time
                              %step (0=no/1=yes) 13) Is the duck wearing a
                              %facemask?
                         
%Pond size in metres (1972 m^2)
Lx = 34.0;
Ly = 58.0;
%time-step length
dt = 1;
%time for duck to recover to full speed
recovery_time = 600;
%Array used to attribute colours to agents based on infected/immune
%properties
colorID = zeros(N);
Depp_sim = 1; %change 0/1 to simulate using ducks with facemasks
%set N initial duck [x,y] position
duck_properties(:,1) = Lx*rand(1,N);
duck_properties(:,2) = Ly*rand(1,N);


%Random number to see if duck will stop/wait leave water
for i = 3:4
    duck_properties(:,i)=rand(1,N);
end

%Set one random duck as infected/immune
b = randi(N);
c = randi(N);
%loop to stop infected duck from being immune duck in initial state 
while c==b    
    c = randi(N);
end
duck_properties(b,6)=1; %infected
duck_properties(c,8)=1; %immuned
%Runs ducks with facemasks sets 50 ducks to have facemasks
if Depp_sim==1
    for i=1:(51)
        duck_properties(i,13) = 1;
    end
end

%Set duck direction & Update duck position
for j=2:TMAX
    %if 10 minutes pass re-randomise the duck stop probability.
    if rem(j,(10*60))==1
        duck_properties(:,3)=rand(1,N);
        duck_properties(:,4)=rand(1,N);
    end
    %loop through agents to update position etc. 
    for i = 1:N
        a=0;
        %if immune colour = dark blue
        if duck_properties(i,8) == 1 
            colorID(i) = 0;        
        %if not infected or immune colour, light blue
        elseif duck_properties(i,6) == 0 && duck_properties(i,8) ~= 1 && duck_properties(i,13)==0
            colorID(i)= 1;
        %if not infected or immune colour and facemask, olive green
        elseif duck_properties(i,6) == 0 && duck_properties(i,8) ~= 1 && duck_properties(i,13)==1
            colorID(i)= 2;
        %If infected colour yellow
        else
            colorID(i)=3;
            %Tracks the time a duck has been infected for (maxes out at 601 seconds)
            if duck_properties(i,7)<= recovery_time
                duck_properties(i,7) = duck_properties(i,7)+1;
            %Sets the immune state of a duck
            else
                if duck_properties(i,8)==1
                    duck_properties(i,7)=0;
                    duck_properties(i,8)=1;
                end 
            end
        end
        %set random direction of the duck
        duck_properties(:,5) = 2*pi*(rand(1,N)-0.5);
        %Sets a fraction of the ducks to remain stationary and tracks the
        %time it has remained stationary
        if duck_properties(i,3)<0.33 && duck_properties(i,9)<10*60
            duck_properties(i,1) = duck_properties(i,1);
            duck_properties(i,2) = duck_properties(i,2);
            duck_properties(i,9) = duck_properties(i,9)+1;
            
        %Setting up boundary conditions. Including the probablilistic
        %nature of grounding a duck and re-entering the pond.
        elseif (duck_properties(i,1)<=0 || duck_properties(i,1)>=Lx || duck_properties(i,2)<=0 || duck_properties(i,2)>=Ly)
            if duck_properties(i,4)<0.1 && duck_properties(i,10)< 60*60
                duck_properties(i,1) = duck_properties(i,1);
                duck_properties(i,2) = duck_properties(i,2);
                duck_properties(i,10) = duck_properties(i,10)+1;
            %if duck ends up outside the box in -x, bounce back right into
            %pond and reset grounded timer
            elseif (duck_properties(i,4)>0.1 || duck_properties(i,10)>3600)  && (duck_properties(i,1)<=0) 
                %if still recovering bounce back half the amount a
                %recovered/normal duck would
                if duck_properties(i,6)==1 && duck_properties(i,7)<=recovery_time
                    duck_properties(i,1) = duck_properties(i,1)+0.15;
                    duck_properties(i,10) = 0;
                else
                    duck_properties(i,1) = duck_properties(i,1)+0.3;
                    duck_properties(i,10) = 0;
                end
            %if duck ends up outside the box in -y, bounce back upwards into
            %pond and reset grounded timer
            elseif (duck_properties(i,4)>0.1 || duck_properties(i,10)>3600)  && (duck_properties(i,2)<=0) 
                if duck_properties(i,6)==1 && duck_properties(i,7)<=recovery_time
                    duck_properties(i,2) = duck_properties(i,2)+0.15;
                    duck_properties(i,10) = 0;
                else
                    duck_properties(i,2) = duck_properties(i,2)+0.3;
                    duck_properties(i,10) = 0;
                end
            %if duck ends up outside the box in +x, bounce back left into
            %pond and reset grounded timer
            elseif (duck_properties(i,4)>0.1 || duck_properties(i,10)>3600)  && (duck_properties(i,1)>=Lx) 
                if duck_properties(i,6)==1 && duck_properties(i,7)<=recovery_time
                    duck_properties(i,1) = duck_properties(i,1)-0.15;
                    duck_properties(i,10) = 0;
                else
                    duck_properties(i,1) = duck_properties(i,1)-0.3;
                    duck_properties(i,10) = 0;
                end
            %if duck ends up outside the box in +y, bounce back downward into
            %pond and reset grounded timer
            elseif (duck_properties(i,4)>0.1 || duck_properties(i,10)>3600)  && (duck_properties(i,2)>=Ly) 
                if duck_properties(i,6)==1 && duck_properties(i,7)<=recovery_time
                    duck_properties(i,2) = duck_properties(i,2)-0.15;
                    duck_properties(i,10) = 0;
                else
                    duck_properties(i,2) = duck_properties(i,2)-0.3;
                    duck_properties(i,10) = 0;
                end
            end        
        else
            %Updates the position of the duck based on whether it has
            %been infected recently
            if duck_properties(i,7)<=recovery_time && duck_properties(i,6)==1
                duck_properties(i,1) = duck_properties(i,1)+0.5*v*cos(duck_properties(i,5));
                duck_properties(i,2) = duck_properties(i,2)+0.5*v*sin(duck_properties(i,5));
                duck_properties(i,9) = 0;
            else
                %normal motion of duck
                duck_properties(i,1) = duck_properties(i,1)+v*cos(duck_properties(i,5));
                duck_properties(i,2) = duck_properties(i,2)+v*sin(duck_properties(i,5));
                duck_properties(i,9) = 0;
            end
        end
    end
    %Formalism for finding distance betwwen ducks
    D = pdist([duck_properties(:,1) duck_properties(:,2)], 'euclidean');
    Z = squareform(D);
    %nested for loop to check if duck is infected/immune. Ducks made
    %infectious/immune cannot retransmit infection/immunity in the timestep
    %it was changed.
    for m=1:100 
        for n=1:100 
            if Z(m,n) < r_transmission %if pairwise distance is smaller than infection radius
                %Checks to determine which way the virus/immunity is spread
                if (duck_properties(m,6)==1 && duck_properties(n,6)~=1 && duck_properties(n,8)~=1) && m ~= n && duck_properties(m,11)==0
                    if Depp_sim==0 || (Depp_sim==1 && (xor(duck_properties(m,13), duck_properties(n,13))==1) && Z(m,n)<=r_transmission/2) || (Depp_sim==1 && duck_properties(m,13)==0 && duck_properties(n,13)==0) || (Depp_sim==1 && duck_properties(m,13)==1 && duck_properties(n,13)==1 && Z(m,n)<=r_transmission/4) 
                        duck_properties(n,6)=1; %set infected
                        duck_properties(n,11)=1; %set infected this time step 
                        %fprintf('Transmission between ducks %d & %d @ time = %d\n', m,n,j)
                    end
                elseif (duck_properties(n,6)==1 && duck_properties(m,6)~=1 && duck_properties(m,8)~=1) && m ~= n && duck_properties(n,11)==0
                    if Depp_sim==0 || (Depp_sim==1 && (xor(duck_properties(m,13), duck_properties(n,13))==1) && Z(m,n)<=r_transmission/2) || (Depp_sim==1 && duck_properties(m,13)==0 && duck_properties(n,13)==0) || (Depp_sim==1 && duck_properties(m,13)==1 && duck_properties(n,13)==1 && Z(m,n)<=r_transmission/4) 
                        duck_properties(m,6)=1; %set infected
                        duck_properties(m,11)=1; %set infected this time step 
                        %fprintf('Transmission between ducks %d & %d @ time = %d\n', m,n,j)
                    end                   
                elseif (duck_properties(m,8)==1 && duck_properties(n,6)==1 && duck_properties(n,7)>=recovery_time) && m ~= n && duck_properties(m,12)==0
                    duck_properties(n,8) = 1; %set immune 
                    duck_properties(n,6) = 0; %sets duck as non-infectious
                    duck_properties(n,12)=1; %set immunised this time step 
                elseif (duck_properties(m,6)==1 && duck_properties(n,8)==1 && duck_properties(m,7)>=recovery_time)  && m ~= n && duck_properties(n,12)==0
                    duck_properties(m,8) = 1; %set immune
                    duck_properties(m,6) = 0; %sets duck as non-infectious
                    duck_properties(m,12)=1; %set immunised this time step 
                end
            end
        end
        %store information about the number of infected/immune/non-infected
        %or non-immune ducks, can output 3/4 graphs dependant on
        %depp_sim=0/1
        if duck_properties(m,6)==1
            plots(sims,1,j) = plots(sims,1,j)+1;
        elseif duck_properties(m,8)==1
            plots(sims,2,j) = plots(sims,2,j)+1;
        elseif duck_properties(m,6)==0 && duck_properties(m,8)==0
            if Depp_sim == 0
                plots(sims,3,j) = plots(sims,3,j)+1;
            elseif Depp_sim==1
                if duck_properties(m,13)==0
                    plots(sims,3,j) = plots(sims,3,j)+1;
                elseif duck_properties(m,13)==1
                    plots(sims,4,j) = plots(sims,4,j)+1;
                end
            end
        end
    end
    %Reset all duck's "infected/immunised this time-step" property (allows a duck to transmit infection disease next time-step)
    duck_properties(:,11)=0;
    duck_properties(:,12)=0;
    %Plotting position of duck after each time-step, take snapshots at
    %t=2/59/600 seconds
    if j==2||j==60||j==recovery_time
        disp(j)
        figure(1)
        scatter(duck_properties(:,1),duck_properties(:,2),[],colorID(:,1))
        %plot circle around patient 0. Rectangular axis scale: circle->elipse.
        p=viscircles([duck_properties(b,1) duck_properties(b,2)],r_transmission,'Color','y');
        p=viscircles([duck_properties(c,1) duck_properties(c,2)],r_transmission,'Color','b');

        box on
        xlim([-2 Lx+2]); ylim([-2 Ly+2]);

        %axis rectangle

        grid on
        %show pond dimensions
        line([0 Lx],[0 0])
        line([0 Lx],[Ly Ly])
        line([0 0],[0 Ly])
        line([Lx Lx],[0 Ly])
        title('Distribution of ducks')
        xlabel('Position in x')
        ylabel('Position in y')
    end
    pause(0.1) 
end
end
%Plot number of infected/immune/etc on a singular axis.
t=linspace(0,TMAX,TMAX);
disp('finished')
av_plots = zeros(4,TMAX);
%calculate the number of ducks in each category at each time step averaged over the number of simulations
for i=1:TMAX
    for j=1:sims
        av_plots(1,i) = av_plots(1,i) + plots(j,1,i);
        av_plots(2,i) = av_plots(2,i) + plots(j,2,i);
        av_plots(3,i) = av_plots(3,i) + plots(j,3,i);
        av_plots(4,i) = av_plots(4,i) + plots(j,4,i);
    end
end
av_plots(:,:) = av_plots(:,:)./n_sims;
%Plot number of infected/immune/etc on a singular axis.
figure(2)
plot(t,av_plots(1,:),'DisplayName', ['Infected'])
hold on
plot(t,av_plots(2,:),'DisplayName', ['Immunised'])
if Depp_sim==0
    plot(t,av_plots(3,:),'DisplayName', ['Non-infected/-immunised'])
end
if Depp_sim==1
    plot(t,av_plots(3,:),'DisplayName', ['Non-infected/-immunised'])
    plot(t,av_plots(4,:),'DisplayName', ['Non-infected/-immunised w/ PPE'])
end

xlim([0 TMAX])
legend('show')

