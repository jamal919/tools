  /*-----------------------------------------------------------------------------
    writing mean tracks time series (track-ref.POC.TP.*.dat file) */

  sprintf(line,"%s.premier_tiers",alti_file);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track data file");
 
  write_list_header(meantrk_file,nXpoints-idum);
  for(n=0,i=0;n<nXpoints;n++) {
      if(projected[n].count==0) continue;
      i++;
      fprintf(meantrk_file,"Pt  : %d\n",i);
      fprintf(meantrk_file,"lon : %lf\n",nominal_track->x[n]);
      fprintf(meantrk_file,"lat : %lf\n",nominal_track->y[n]);
      fprintf(meantrk_file,"Mes : %d\n",(int)(projected[n].count/3.0));
         
      for(k=0;k<(int)(projected[n].count/3.0);k++) {
          fprintf(meantrk_file,"%12.6f",projected[n].data[k].time);
          fprintf(meantrk_file," %7.4f",projected[n].data[k].values[28]-projected[n].data[k].values[1]-projected[n].data[k].values[10]-projected[n].data[k].values[11]);
          fprintf(meantrk_file,"\n");
        }
      fprintf(meantrk_file,"#-----------------------------------------------------\n");
    }   /* loop on nXpoints */
  fflush(meantrk_file);
  fclose(meantrk_file);
   /*-----------------------------------------------------------------------------
    writing mean tracks time series (track-ref.POC.TP.*.dat file) */

  sprintf(line,"%s.second_tiers",alti_file);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track data file");
 
  write_list_header(meantrk_file,nXpoints-idum);
  for(n=0,i=0;n<nXpoints;n++) {
      if(projected[n].count==0) continue;
      i++;
      fprintf(meantrk_file,"Pt  : %d\n",i);
      fprintf(meantrk_file,"lon : %lf\n",nominal_track->x[n]);
      fprintf(meantrk_file,"lat : %lf\n",nominal_track->y[n]);
      int toto=(int)floor(projected[n].count-projected[n].count/3.0*2.0);
      fprintf(meantrk_file,"Mes : %d\n",toto);
         
      for(k=0;k<toto;k++) {
          fprintf(meantrk_file,"%12.6f",projected[n].data[k+toto].time);
          fprintf(meantrk_file," %7.4f",projected[n].data[k+toto].values[28]-projected[n].data[k+toto].values[1]-projected[n].data[k+toto].values[10]-projected[n].data[k+toto].values[11]);
          fprintf(meantrk_file,"\n");
        }
      fprintf(meantrk_file,"#-----------------------------------------------------\n");
    }   /* loop on nXpoints */
  fflush(meantrk_file);
  fclose(meantrk_file);

 /*-----------------------------------------------------------------------------
    writing mean tracks time series (track-ref.POC.TP.*.dat file) */

  sprintf(line,"%s.dernier_tiers",alti_file);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track data file");
 
  write_list_header(meantrk_file,nXpoints-idum);
  for(n=0,i=0;n<nXpoints;n++) {
      if(projected[n].count==0) continue;
      i++;
      fprintf(meantrk_file,"Pt  : %d\n",i);
      fprintf(meantrk_file,"lon : %lf\n",nominal_track->x[n]);
      fprintf(meantrk_file,"lat : %lf\n",nominal_track->y[n]);
      int toto=(int)floor(projected[n].count-projected[n].count/3.0*2.0);
      fprintf(meantrk_file,"Mes : %d\n",toto);
         
      for(k=0;k<toto;k++) {
          fprintf(meantrk_file,"%12.6f",projected[n].data[k+toto+toto].time);
          fprintf(meantrk_file," %7.4f",projected[n].data[k+toto+toto].values[28]-projected[n].data[k+toto+toto].values[1]-projected[n].data[k+toto+toto].values[10]-projected[n].data[k+toto+toto].values[11]);
          fprintf(meantrk_file,"\n");
        }
      fprintf(meantrk_file,"#-----------------------------------------------------\n");
    }   /* loop on nXpoints */
  fflush(meantrk_file);
  fclose(meantrk_file);

  /*-----------------------------------------------------------------------------
    writing mean tracks time series (track-ref.POC.TP.*.dat file) */

  sprintf(line,"%s.FES2004_analyse_list",alti_file);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track data file");
 
  write_list_header(meantrk_file,nXpoints-idum);
  for(n=0,i=0;n<nXpoints;n++) {
      if(projected[n].count==0) continue;
      i++;
      fprintf(meantrk_file,"Pt  : %d\n",i);
      fprintf(meantrk_file,"lon : %lf\n",nominal_track->x[n]);
      fprintf(meantrk_file,"lat : %lf\n",nominal_track->y[n]);
      fprintf(meantrk_file,"Mes : %d\n",projected[n].count);
         
      for(k=0;k<projected[n].count;k++) {
          fprintf(meantrk_file,"%12.6f",projected[n].data[k].time);
          // fprintf(meantrk_file," %7.4f",projected[n].data[k].values[28]-projected[n].data[k].values[1]-projected[n].data[k].values[10]-projected[n].data[k].values[11]);
          fprintf(meantrk_file,"       %7.4f",projected[n].data[k].values[18]);
          fprintf(meantrk_file,"\n");
        }
      fprintf(meantrk_file,"#-----------------------------------------------------\n");
    }   /* loop on nXpoints */
  fflush(meantrk_file);
  fclose(meantrk_file);
  /*-----------------------------------------------------------------------------
    writing mean tracks time series (track-ref.POC.TP.*.dat file) */

  sprintf(line,"%s.LOAD_analyse_list",alti_file);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track data file");
 
  write_list_header(meantrk_file,nXpoints-idum);
  for(n=0,i=0;n<nXpoints;n++) {
      if(projected[n].count==0) continue;
      i++;
      fprintf(meantrk_file,"Pt  : %d\n",i);
      fprintf(meantrk_file,"lon : %lf\n",nominal_track->x[n]);
      fprintf(meantrk_file,"lat : %lf\n",nominal_track->y[n]);
      fprintf(meantrk_file,"Mes : %d\n",projected[n].count);
         
      for(k=0;k<projected[n].count;k++) {
          fprintf(meantrk_file,"%12.6f",projected[n].data[k].time);
          // fprintf(meantrk_file," %7.4f",projected[n].data[k].values[28]-projected[n].data[k].values[1]-projected[n].data[k].values[10]-projected[n].data[k].values[11]);
          fprintf(meantrk_file,"       %7.4f",projected[n].data[k].values[10]);
          fprintf(meantrk_file,"\n");
        }
      fprintf(meantrk_file,"#-----------------------------------------------------\n");
    }   /* loop on nXpoints */
  fflush(meantrk_file);
  fclose(meantrk_file);
  /*-----------------------------------------------------------------------------
    writing mean tracks time series (track-ref.POC.TP.*.dat file) */

  sprintf(line,"%s.geocentri_analyse_list",alti_file);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track data file");
 
  write_list_header(meantrk_file,nXpoints-idum);
  for(n=0,i=0;n<nXpoints;n++) {
      if(projected[n].count==0) continue;
      i++;
      fprintf(meantrk_file,"Pt  : %d\n",i);
      fprintf(meantrk_file,"lon : %lf\n",nominal_track->x[n]);
      fprintf(meantrk_file,"lat : %lf\n",nominal_track->y[n]);
      fprintf(meantrk_file,"Mes : %d\n",projected[n].count);
         
      for(k=0;k<projected[n].count;k++) {
          fprintf(meantrk_file,"%12.6f",projected[n].data[k].time);
          // fprintf(meantrk_file," %7.4f",projected[n].data[k].values[28]-projected[n].data[k].values[1]-projected[n].data[k].values[10]-projected[n].data[k].values[11]);
          fprintf(meantrk_file,"       %7.4f",projected[n].data[k].values[28]);
          fprintf(meantrk_file,"\n");
        }
      fprintf(meantrk_file,"#-----------------------------------------------------\n");
    }   /* loop on nXpoints */
  fflush(meantrk_file);
  fclose(meantrk_file);
  /*-----------------------------------------------------------------------------
    writing mean tracks time series (track-ref.POC.TP.*.dat file) */

  sprintf(line,"%s.mog2d_analyse_list",alti_file);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track data file");
 
  write_list_header(meantrk_file,nXpoints-idum);
  for(n=0,i=0;n<nXpoints;n++) {
      if(projected[n].count==0) continue;
      i++;
      fprintf(meantrk_file,"Pt  : %d\n",i);
      fprintf(meantrk_file,"lon : %lf\n",nominal_track->x[n]);
      fprintf(meantrk_file,"lat : %lf\n",nominal_track->y[n]);
      fprintf(meantrk_file,"Mes : %d\n",projected[n].count);
         
      for(k=0;k<projected[n].count;k++) {
          fprintf(meantrk_file,"%12.6f",projected[n].data[k].time);
          // fprintf(meantrk_file," %7.4f",projected[n].data[k].values[28]-projected[n].data[k].values[1]-projected[n].data[k].values[10]-projected[n].data[k].values[11]);
          fprintf(meantrk_file,"       %7.4f",projected[n].data[k].values[1]);
          fprintf(meantrk_file,"\n");
        }
      fprintf(meantrk_file,"#-----------------------------------------------------\n");
    }   /* loop on nXpoints */
  fflush(meantrk_file);
  fclose(meantrk_file);
  /*-----------------------------------------------------------------------------
    writing mean tracks time series (track-ref.POC.TP.*.dat file) */

  sprintf(line,"%s.solid_analyse_list",alti_file);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track data file");
 
  write_list_header(meantrk_file,nXpoints-idum);
  for(n=0,i=0;n<nXpoints;n++) {
      if(projected[n].count==0) continue;
      i++;
      fprintf(meantrk_file,"Pt  : %d\n",i);
      fprintf(meantrk_file,"lon : %lf\n",nominal_track->x[n]);
      fprintf(meantrk_file,"lat : %lf\n",nominal_track->y[n]);
      fprintf(meantrk_file,"Mes : %d\n",projected[n].count);
         
      for(k=0;k<projected[n].count;k++) {
          fprintf(meantrk_file,"%12.6f",projected[n].data[k].time);
          // fprintf(meantrk_file," %7.4f",projected[n].data[k].values[28]-projected[n].data[k].values[1]-projected[n].data[k].values[10]-projected[n].data[k].values[11]);
          fprintf(meantrk_file,"       %7.4f",projected[n].data[k].values[11]);
          fprintf(meantrk_file,"\n");
        }
      fprintf(meantrk_file,"#-----------------------------------------------------\n");
    }   /* loop on nXpoints */
  fflush(meantrk_file);
  fclose(meantrk_file);
