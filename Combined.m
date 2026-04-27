%Combined Program

%Welcome Prompts
disp(sprintf ('Welcome to ReactM'));
disp(sprintf ('Available Models Listed'));

%Reactor Models
disp( sprintf('1. Batch Reactor'));
disp( sprintf('2. Continuous Stirred Tank Reactor'));
disp( sprintf('3. Plug Flow Reactor'));
disp( sprintf('4. Packed Bed Reactor'));
disp( sprintf('5. Rigorous CSTR'));
disp( sprintf('6. Rigorous PFR'));
disp( sprintf(' '));
             
prompt = ' Select model to run   ';
simulation = input(prompt);

if (simulation == 1)

    % Batch
disp('Running Batch Reactor');
open('Batch.m'); % Execuion 

elseif (simulation == 2)

    % Continuous Stirred Tank Reactor
disp('Running Continuous Stirred Tank Reactor');
open("ContinuousStirredTankReactor.m"); % Execution


elseif (simulation == 3)

    % Plug Flow Reactor
disp('Running Plug Flow Reactor');
open("PlugFlowReactor.m"); % Execution

elseif (simulation == 4)

    % Packed Bed Reactor
disp('Packed Bed Reactor');
open("PackedBedReactor.m"); % Execution

elseif (simulation == 5)

    % Rigorous CSTR
disp('Rigorous CSTR');
open('CSTRPRG.m'); % Execution

elseif (simulation == 6)

    % Rigorous PFR
disp('Rigorous PFR');
open("PFRPRG.m"); % Execution

end

