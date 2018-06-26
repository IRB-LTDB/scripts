%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example query to the LTDB database                                      %
% Notes:                                                                  %
% - Requires the mysql-connector-java-5.x.xx-bin.jar in the same folder   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;
javaaddpath mysql-connector-java-5.1.45-bin.jar

% SETTINGS
SERVER = 'localhost'
DB_NAME = 'ltdb';
USER = 'root';
PASSWORD = '';
DBMS_TYPE = 'MySQL';

% Create a connection to the database
conn = database(DB_NAME, USER, PASSWORD, 'Vendor', DBMS_TYPE, 'Server', SERVER);

% Retrieves the ID of the tracks of Neutrophils only
strQuery = ['SELECT t.id_track FROM tracks t, cell_types ct ', ...
            'WHERE t.id_cell_type = ct.id_cell_type ', ...
            'AND ct.descr = ''Neutrophils'''];
        
% Executes the query
data = select(conn, strQuery);

close(conn);
        
