
#include <qapplication.h>
#include <qplatinumstyle.h>

#include "MeshValmetControls.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
using namespace std;

void print_usage(FILE *out)
{
  fprintf(out,"Usage: MeshValmet [options]\n");
  fprintf(out,"options:");
  fprintf(out,"\n");
  fprintf(out,"  -p\tbatchFile\n");
  fprintf(out,"  -s\tsamplingStep\t(Default value = 0.5)\n");
  fprintf(out,"  -f\tMinSamplingFreq\t(Default value = 2)\n");
  fprintf(out,"  -b\tNumOfBins\t(Default value = 256)\n");
  fprintf(out,"  -type\t computeType\t(Use 0 to compute the error of A <->B,\n");
  fprintf(out,"  \t\t\tUse 1 to compute the error of A ->B,\n");
  fprintf(out,"  \t\t\tUse 2 to compute the error of B ->A,\n");
  fprintf(out,"  \t\t\tDefault value = 0)\n");
  fprintf(out,"  -absolute\t\t(Default FALSE)\n");
  fprintf(out,"  -in1\tinFile1\n");
  fprintf(out,"  -in2\tinFile2\n");
  fprintf(out,"  -o\toutFile\n");   
  fprintf(out,"\n");
  fprintf(out,"Batch file format:\n");
  fprintf(out,"  numOfFiles\n");
  fprintf(out,"  inFile11\tinFile12\toutFile1\n");
  fprintf(out,"  inFile21\tinFile22\toutFile2\n");
  fprintf(out,"  ......\t......\t......\n");
}

/* Output: Max, Min, Vertex Mean, Vertex Sigma, Face Mean, Face RMS,
 * Hausdorff, MSD, MAD, Median, 68percentile1, 68percentile2,
 * 95perentile1, 95percentile2, 75percentile1, 75percentile2 dist volume
 */
int compare(char* in1, char* in2, char* out, float samStep, int frequency, int type, int bins, bool absolute)
{
  //cout<<samStep<<":"<<frequency<<":"<<type<<":"<<bins<<":"<<absolute<<endl;
  
  struct outbuf *log;
  
  struct model_error model1,model2;
  struct args pargs;
  struct dist_surf_surf_stats stats;
  struct dist_surf_surf_stats stats_rev;
  double abs_sampling_step, abs_sampling_dens;
  
  int *histogram1;
  int *histogram2;
  int *histogram3;
  int *hist;
  
  double hausdorff_error;
  double abs_sum_error, sum_sqr_error, mad_error, msd_error, dist_volume;
  //mean_error is vertex mean; sigma_error is vertex sigma(rms).
  double mean_error,sum_error,sum_square_deviations, sigma_error, error_max, error_min;
  //face_mean_error is mean of face error; face_rms_eroor is face sigma(rms).
  double face_mean_error, face_rms_error;
  //Median; 95percentile; 68precentile
  double median, percentile95_1, percentile95_2,percentile68_1, percentile68_2, percentile75_1, percentile75_2;    

  memset(&model1,0,sizeof(model1));
  memset(&model2,0,sizeof(model2));
  log = NULL;

  //mesh_run: mesh_run(&pargs, &model1, &model2, log, &pr, &stats, &stats_rev, &abs_sampling_step, &abs_sampling_dens);
  pargs.m1_fname =  in1;
  pargs.m2_fname =  in2;
  pargs.sampling_step = samStep*0.01;
  pargs.min_sample_freq = frequency;
  pargs.do_wlog = false;
  log = outbuf_new(stdio_puts,stdout);
  pargs.verb_analysis = false;
  pargs.do_texture = false;
  if(type == 0)
    pargs.do_symmetric = 3;
  else if(type == 1)
    pargs.do_symmetric = 1;
  else if(type == 2)
    pargs.do_symmetric = 2;  
  pargs.no_gui = 1;
  
  mesh_run(&pargs, &model1, &model2, log, NULL, &stats, &stats_rev, &abs_sampling_step, &abs_sampling_dens);

  double *serror;
  int i,n;
  abs_sum_error = 0;
  sum_sqr_error = 0;
  sum_error = 0;
  error_max = -100000;
  error_min = 100000;
  sum_square_deviations = 0;  
  
  double drange,off;
  int bin_idx;
  int len = bins;
  
  // Matthieu FIX: Don't delete histograms here
  /*delete [] histogram1;
  delete [] histogram2;
  delete [] hist;*/

  histogram1 = new int[len];
  histogram2 = new int[len];
  histogram3 = new int[len];
  hist = new int[len];
  memset(histogram1, 0, sizeof(*histogram1)*len);
  memset(histogram2, 0, sizeof(*histogram2)*len);
  memset(histogram3, 0, sizeof(*histogram3)*len);
  memset(hist, 0, sizeof(*hist)*len);

  //do the stat
  if(absolute)
  {
    if(type == 0)
    {
      int n1, n2;
      n1 = model1.n_samples;
      serror = model1.fe[0].serror;
      // We calculate the mean and the SD (sigma) error
      for (i=0; i<n1; i++) 
      { 
        sum_error += fabs(serror[i]);
        sum_sqr_error += serror[i]*serror[i];
        if(error_max < fabs(serror[i]))
          error_max = fabs(serror[i]);
        if(error_min > fabs(serror[i]))
          error_min = fabs(serror[i]);
      }
      n2 = model2.n_samples;
      serror = model2.fe[0].serror;
      
      for (i=0; i<n2; i++) 
      { 
        sum_error += fabs(serror[i]);
        sum_sqr_error += serror[i]*serror[i];
        if(error_max < fabs(serror[i]))
          error_max = fabs(serror[i]);
        if(error_min > fabs(serror[i]))
          error_min = fabs(serror[i]);
      }
      n = n1 + n2;
      mean_error = sum_error / n;
      mad_error = sum_error / n;
      msd_error = sum_sqr_error /n;  
          
      serror = model1.fe[0].serror;
      for (i=0; i<n1; i++) 
            sum_square_deviations += (fabs(serror[i])-mean_error)*(fabs(serror[i])-mean_error);
      serror = model2.fe[0].serror;
      for (i=0; i<n2; i++) 
          sum_square_deviations += (fabs(serror[i])-mean_error)*(fabs(serror[i])-mean_error);
    
      sigma_error = sqrt(sum_square_deviations/(n-1));
        
      face_mean_error = (stats.abs_mean_tot+stats_rev.abs_mean_tot)/(stats.st_m1_area+stats_rev.st_m1_area);
      dist_volume = (stats.abs_mean_tot+stats_rev.abs_mean_tot)/2.0;
      face_rms_error = sqrt((stats.abs_rms_tot+stats_rev.abs_rms_tot)/(stats.st_m1_area+stats_rev.st_m1_area));
    }
    else if(type == 1)
    {
      n = model1.n_samples;
      cout<<n<<endl;
      serror = model1.fe[0].serror;
      // We calculate the mean and the SD (sigma) error
      for (i=0; i<n; i++) 
      { 
        sum_error += fabs(serror[i]);
        sum_sqr_error += serror[i]*serror[i];
        if(error_max < fabs(serror[i]))
          error_max = fabs(serror[i]);
        if(error_min > fabs(serror[i]))
          error_min = fabs(serror[i]);
      }
      mean_error = sum_error / n;
      mad_error = sum_error / n;
      msd_error = sum_sqr_error /n;

      for (i=0; i<n; i++) 
        sum_square_deviations += (fabs(serror[i])-mean_error)*(fabs(serror[i])-mean_error);
    
      sigma_error = sqrt(sum_square_deviations/(n-1));
        
      face_mean_error = stats.abs_mean_tot / stats.st_m1_area;
      dist_volume = stats.abs_mean_tot;
      face_rms_error = stats.abs_rms_dist;
    }
    else if(type == 2)
    {
      n = model2.n_samples;
      serror = model2.fe[0].serror;
      // We calculate the mean and the SD (sigma) error
      for (i=0; i<n; i++) 
      { 
        sum_error += fabs(serror[i]);
        sum_sqr_error += serror[i]*serror[i];
        if(error_max < fabs(serror[i]))
          error_max = fabs(serror[i]);
        if(error_min > fabs(serror[i]))
          error_min = fabs(serror[i]);
      }
      mean_error = sum_error / n;
      mad_error = sum_error / n;
      msd_error = sum_sqr_error /n;
          
      for (i=0; i<n; i++) 
           sum_square_deviations += (fabs(serror[i])-mean_error)*(fabs(serror[i])-mean_error);
    
      sigma_error = sqrt(sum_square_deviations/(n-1));
        
      face_mean_error = stats_rev.abs_mean_tot / stats_rev.st_m1_area;
      dist_volume = stats_rev.abs_mean_tot;
      face_rms_error = stats_rev.abs_rms_dist;
    }
    
    hausdorff_error = error_max;  
    
    //compute histogram
    if(type == 1 || type == 0)
    {
        n = model1.n_samples;
        drange = error_max-error_min;
        off = error_min;
        serror = model1.fe[0].serror;
        for (i=0; i<n; i++) 
        {
           bin_idx = (int) ((fabs(serror[i])-off)/drange*len);
            if (bin_idx >= len) bin_idx = len-1;
            histogram1[bin_idx]++;
        }
    }
    if(type == 2 || type == 0)
    {
        n = model1.n_samples;
        drange = error_max-error_min;
        off = error_min;
        serror = model2.fe[0].serror;
        for (i=0; i<n; i++) 
        {
           bin_idx = (int) ((fabs(serror[i])-off)/drange*len);
          if (bin_idx >= len) bin_idx = len-1;
            histogram2[bin_idx]++;
        }  
    }
    if(type == 0)
    {
      for(i=0; i<len; i++)
      hist[i] = histogram1[i] + histogram2[i];  
    }
  }
  else
  {
    if(type == 0)
    {
      int n1, n2;
      n1 = model1.n_samples;
      serror = model1.fe[0].serror;
      // We calculate the mean and the SD (sigma) error
      for (i=0; i<n1; i++) 
      { 
        sum_error += serror[i];
        abs_sum_error += fabs(serror[i]);
        sum_sqr_error += serror[i]*serror[i];
        if(error_max < serror[i])
          error_max = serror[i];
        if(error_min > serror[i])
          error_min = serror[i];
      }
      n2 = model2.n_samples;
      serror = model2.fe[0].serror;
      
      for (i=0; i<n2; i++) 
      { 
        sum_error += serror[i];
        abs_sum_error += fabs(serror[i]);
        sum_sqr_error += serror[i]*serror[i];
        if(error_max < serror[i])
          error_max = serror[i];
        if(error_min > serror[i])
          error_min = serror[i];
      }
      n = n1 + n2;
      mean_error = sum_error / n;
      mad_error = abs_sum_error / n;
      msd_error = sum_sqr_error /n;  
          
      serror = model1.fe[0].serror;
      for (i=0; i<n1; i++) 
            sum_square_deviations += (serror[i]-mean_error)*(serror[i]-mean_error);
      serror = model2.fe[0].serror;
      for (i=0; i<n2; i++) 
          sum_square_deviations += (serror[i]-mean_error)*(serror[i]-mean_error);
    
      sigma_error = sqrt(sum_square_deviations/(n-1));

      face_mean_error = (stats.mean_tot+stats_rev.mean_tot)/(stats.st_m1_area+stats_rev.st_m1_area);
      dist_volume = (stats.abs_mean_tot+stats_rev.abs_mean_tot)/2.0;
      face_rms_error = sqrt((stats.rms_tot+stats_rev.rms_tot)/(stats.st_m1_area+stats_rev.st_m1_area));
    }
    else if(type == 1)
    {
      n = model1.n_samples;
      serror = model1.fe[0].serror;
      // We calculate the mean and the SD (sigma) error
      for (i=0; i<n; i++) 
      { 
        sum_error += serror[i];
        abs_sum_error += fabs(serror[i]);
        sum_sqr_error += serror[i]*serror[i];
        if(error_max < serror[i])
          error_max = serror[i];
        if(error_min > serror[i])
          error_min = serror[i];
      }
      mean_error = sum_error / n;
      mad_error = abs_sum_error / n;
      msd_error = sum_sqr_error /n;

      for (i=0; i<n; i++) 
        sum_square_deviations += (serror[i]-mean_error)*(serror[i]-mean_error);
    
      sigma_error = sqrt(sum_square_deviations/(n-1));

      face_mean_error = stats.mean_dist;
      dist_volume = stats.abs_mean_tot;
      face_rms_error = stats.rms_dist;
    }
    else if(type == 2)
    {
      n = model2.n_samples;
      serror = model2.fe[0].serror;
      // We calculate the mean and the SD (sigma) error
      for (i=0; i<n; i++) 
      { 
        sum_error += serror[i];
        abs_sum_error += fabs(serror[i]);
        sum_sqr_error += serror[i]*serror[i];
        if(error_max < serror[i])
          error_max = serror[i];
        if(error_min > serror[i])
          error_min = serror[i];
      }
      mean_error = sum_error / n;
      mad_error = abs_sum_error / n;
      msd_error = sum_sqr_error /n;
          
      for (i=0; i<n; i++) 
           sum_square_deviations += (serror[i]-mean_error)*(serror[i]-mean_error);
    
      sigma_error = sqrt(sum_square_deviations/(n-1));

      face_mean_error = stats_rev.mean_dist;
      dist_volume = stats_rev.abs_mean_tot;
      face_rms_error = stats_rev.rms_dist;
    }
    
    hausdorff_error = fabs(error_max)>fabs(error_min)?fabs(error_max):fabs(error_min);  
    
    //compute histogram
    if(type == 1 || type == 0)
    {
        n = model1.n_samples;
        drange = model1.max_error-model1.min_error;
        off = model1.min_error;
        serror = model1.fe[0].serror;
        for (i=0; i<n; i++) 
        {
           bin_idx = (int) ((serror[i]-off)/drange*len);
            if (bin_idx >= len) bin_idx = len-1;
            histogram1[bin_idx]++;
        }
    }
    if(type == 2 || type == 0)
    {
        n = model2.n_samples;
        drange = model2.max_error-model2.min_error;
        off = model2.min_error;
        serror = model2.fe[0].serror;
        for (i=0; i<n; i++) 
        {
           bin_idx = (int) ((serror[i]-off)/drange*len);
          if (bin_idx >= len) bin_idx = len-1;
            histogram2[bin_idx]++;
        }  
    }
    /*
    // cubic root histogram
    n = model1.n_samples;
    
    drange = exp(log(model1.max_error)/3)-exp(log(model1.min_error)/3);
    off = exp(log(model1.min_error)/3);
    serror = exp(log(model1.fe[0].serror)/3);
    for (i=0; i<n; i++) 
    {
       bin_idx = (int) ((exp(log(serror[i])/3)-off)/drange*len);
       if (bin_idx >= len) bin_idx = len-1;
         histogram3[bin_idx]++;
    }
    // end cubic root histogram
    */
    
    if(type == 0)
    {
      for(i=0; i<len; i++)
      {
        hist[i] = histogram1[i] + histogram2[i];  
        //cout<<hist[i]<<"\t";
      }
    }  
  }
  
  //cout<<error_max<<":"<<error_min<<":"<<msd_error<<":"<<mad_error<<":"<<mean_error<<":"<<dist_volume<<endl;
  
  //Calculate the median
  int sum_cnt=0;
  int* curr_hist;
  double half_sample_no;
  
  if(type == 1)
  {
    half_sample_no =  model1.n_samples / 2.0 ;
    curr_hist = histogram1;
  }
  else if(type == 2)
  {
    half_sample_no =  model2.n_samples / 2.0 ;
    curr_hist = histogram2;
  }
  else if(type == 0)
  {
    half_sample_no =  (model1.n_samples + model2.n_samples) / 2.0 ;
    curr_hist = hist;        
  }
  
  for(i = 0; i < len; i++)
  {
      if(sum_cnt >=  half_sample_no )
        break;
      else
        sum_cnt +=  curr_hist[i];
  }
  //cout<<half_sample_no<<":"<<len<<":"<<i<<"\n";

  median = error_min + (error_max - error_min)*(i-1)/(double)len;
  median = median - (error_max-error_min)*(sum_cnt-half_sample_no)/(double)(len*curr_hist[i-1]);
  //cout<<error_min<<":"<<error_max<<":"<<i<<":"<<half_sample_no<<endl;
  //cout<<"meadian:"<<median<<endl;

  //compute percentile
  double percent95= 0.95;
  double percent68 = 0.68;
  double percent75 = 0.75;
  percentile95_1 = 0;
  percentile95_2 = 0;
  percentile68_1 = 0;
  percentile68_2 = 0;
  percentile75_1 = 0;
  percentile75_2 = 0;

  int curr_samples;
  
  if(type == 1)
  {    
    curr_samples = model1.n_samples;
    curr_hist = histogram1;    
  }
  else if(type == 2)
  {
    curr_samples = model2.n_samples;
    curr_hist = histogram2;  
  }
  else if(type == 0)
  {
    curr_samples = model1.n_samples + model2.n_samples;
    curr_hist = hist;
  }
  
  // Previous computation was wrong
  // 95% percentile should compute 5%...95%
  // Official def: The 95th percentile is the smallest number that is 
  //               greater that 95% of the numbers in a given set. 
  // the previous code computed the 95% percentile as the range around the median that covers 95%
  //  I don't really know what that is? but it was way harder to compute...
   
  sum_cnt = 0;
  for(i = 0; i < len; i++)
  {
    sum_cnt +=  curr_hist[i];
    if(sum_cnt >= (double) curr_samples * (1.0 - percent95)  ) {
      percentile95_1 = error_min + (error_max-error_min)*i/(double)len;
      break;
    }     
  }
  sum_cnt = 0;
  for(i = 0; i < len; i++)
  {
    sum_cnt +=  curr_hist[i];
    if(sum_cnt >= (double) curr_samples * percent95  ) {
      percentile95_2 = error_min + (error_max-error_min)*i/(double)len;
      break;
    }     
  }

  sum_cnt = 0;
  for(i = 0; i < len; i++)
  {
    sum_cnt +=  curr_hist[i];
    if(sum_cnt >= (double) curr_samples * (1.0 - percent68)  ) {
      percentile68_1 = error_min + (error_max-error_min)*i/(double)len;
      break;
    }     
  }
  sum_cnt = 0;
  for(i = 0; i < len; i++)
  {
    sum_cnt +=  curr_hist[i];
    if(sum_cnt >= (double) curr_samples * percent68  ) {
      percentile68_2 = error_min + (error_max-error_min)*i/(double)len;
      break;
    }     
  }

  sum_cnt = 0;
  for(i = 0; i < len; i++)
  {
    sum_cnt +=  curr_hist[i];
    if(sum_cnt >= (double) curr_samples * (1.0 - percent75)  ) {
      percentile75_1 = error_min + (error_max-error_min)*i/(double)len;
      break;
    }     
  }
  sum_cnt = 0;
  for(i = 0; i < len; i++)
  {
    sum_cnt +=  curr_hist[i];
    if(sum_cnt >= (double) curr_samples * percent75  ) {
      percentile75_2 = error_min + (error_max-error_min)*i/(double)len;
      break;
    }     
  }

  //output stat
  char* errorfile = new char[strlen(out)+6];
  char* histfile = new char[strlen(out)+6];
  char* statfile = new char[strlen(out)+6];
  strcpy(errorfile,out);
  strcat(errorfile,".errs");
  strcpy(histfile,out);
  strcat(histfile,".hist");
  strcpy(statfile,out);
  strcat(statfile,".stat");

  FILE* fp;
  fp = fopen(errorfile,"a");
  if(type == 1||type == 0)
  {
    fprintf(fp,"%s->%s\n",in1,in2);
    n = model1.n_samples;
    fprintf(fp,"%d\n",n);
    serror = model1.fe[0].serror;
    
    for (i=0; i<n; i++) 
      fprintf(fp,"%f\t",serror[i]);
    fprintf(fp,"\n\n");
  }
  else if(type == 2||type == 0)
  {
    fprintf(fp,"%s->%s\n",in2,in1);
    n = model2.n_samples;
    fprintf(fp,"%d\n",n);
    serror = model2.fe[0].serror;
    
    for (i=0; i<n; i++) 
      fprintf(fp,"%f\t",serror[i]);
    fprintf(fp,"\n\n");
  }
  fclose(fp);
  
  fp = fopen(histfile,"a");
  if(type == 1||type == 0)
  {
    fprintf(fp,"%s->%s\n",in1,in2);
    fprintf(fp,"Min:\t%f\n",error_min);
    fprintf(fp,"Max:\t%f\n",error_max);
    fprintf(fp,"Bins:\t%d\n",bins);
    for(i = 0; i < len; i++)
      fprintf(fp,"%d\n",histogram1[i]);
    fprintf(fp,"\n\n");
  }
  else if(type == 2||type == 0)
  {
    fprintf(fp,"%s->%s\n",in2,in1);
    fprintf(fp,"Min:\t%f\n",error_min);
    fprintf(fp,"Max:\t%f\n",error_max);
    fprintf(fp,"Bins:\t%d\n",bins);
    for(i = 0; i < len; i++)
      fprintf(fp,"%d\t",histogram2[i]);
    fprintf(fp,"\n\n");
  }  
  /*
  //output cubic root histogram
  fprintf(fp,"cubic root\t%s->%s\n",in1,in2);
  fprintf(fp,"Min:\t%f\n",exp(log(error_min)/3));
  fprintf(fp,"Max:\t%f\n",exp(log(error_max)/3));
  fprintf(fp,"Bins:\t%d\n",bins);
  for(i = 0; i < len; i++)
    fprintf(fp,"%d\t",histogram3[i]);
  fprintf(fp,"\n\n");
  //end output cubic root histogram
  */
  fclose(fp);
  
  fp = fopen(statfile,"a");
  //fprintf(fp,"Min\tMax\tVertexMean\tVertexSigma\tFaceMean\tFaceRMS\tHausdorff\tMSD\tMAD\tMedian\tDistVolume\t68percentile1\t68percentile2\t95percentile1\t95percentile2\t75percentile1\t75percentile2\n");
  fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",error_min, error_max, mean_error, sigma_error, face_mean_error, face_rms_error, hausdorff_error, msd_error, mad_error, median, dist_volume, percentile68_1, percentile68_2, percentile95_1,percentile95_2,percentile75_1,percentile75_2);
  //cout<<error_min<<":"<<error_max<<":"<<mean_error<<":"<<sigma_error<<":"<<face_mean_error<<":"<<face_rms_error<<":"<<hausdorff_error<<":"<<mad_error<<":"<<mad_error<<":"<<median<<":"<<dist_volume<<":"<<percentile68_1<<":"<<percentile68_2<<":"<<percentile95_1<<":"<<percentile95_2<<endl;
  fclose(fp);
  delete [] errorfile;
  delete [] histfile;
  delete [] statfile;
  
  delete [] histogram1;
  delete [] histogram2;
  delete [] histogram3;
  
  return 0;
}

int main( int argc, char* argv[] ) 
{
  if(argc <= 1)
  {
    QApplication myApp( argc, argv );
  
    MeshValmetControls m_GUI( 0, 0 );
    myApp.setMainWidget(&m_GUI);
  
    myApp.setStyle( new QPlatinumStyle() );
  
    m_GUI.show();
    myApp.exec();  
  }
  else
  {
    int i = 1, bins = 256, type = 0;
    float ss = 0.5, freq = 2.0;
    char* inFile1 = "in1.byu";
    char* inFile2 = "in2.byu";    
    char* outFile = "out";
    char* batchFile = "batch";
    bool batch = false, absolute = false;
    
    while (i < argc) 
    {
          if (argv[i][0] == '-') { /* Option */
              if (strcmp(argv[i],"-h") == 0) { /* help */
              print_usage(stdout);
              exit(0);
              }
              else if (strcmp(argv[i],"-p") == 0) { /* batchFile */
                i++;
                batchFile = argv[i];
                batch = true;
              }
              else if (strcmp(argv[i],"-in1") == 0) { /* inFile */
                i++;
                inFile1 = argv[i];
                batch = false;
              }
              else if (strcmp(argv[i],"-in2") == 0) { /* inFile */
                i++;
                inFile2 = argv[i];
                batch = false;
              }
              else if (strcmp(argv[i],"-o") == 0) { /* outFile */
                i++;
                outFile = argv[i];
                batch = false;
              }
              else if (strcmp(argv[i],"-s") == 0) { 
                i++;
                ss = atof(argv[i]);
              }
              else if (strcmp(argv[i],"-f") == 0) { 
                i++;
                freq = atoi(argv[i]);
              }
              else if (strcmp(argv[i],"-b") == 0) { 
                i++;
                bins = atoi(argv[i]);
              }
              else if (strcmp(argv[i],"-type") == 0) { 
                i++;
                type = atoi(argv[i]);
              }
              else if (strcmp(argv[i],"-absolute") == 0) { 
                absolute = true;
              }
              else { /* unrecognized option */
                fprintf(stdout,  "ERROR: unknown option in command line, use -h for help\n");
                exit(1);
              }
          }
          else
          {
            fprintf(stdout,  "ERROR: unknown option in command line, use -h for help\n");
            exit(1);  
          }
          i++;
    }
    //cout<<ss<<":"<<freq<<":"<<bins<<":"<<type<<":"<<absolute<<":"<<inFile1<<":"<<inFile2<<":"<<outFile<<endl;
    
    // Remeber to create space for the '\0'.
    char* errorfile = new char[strlen(outFile)+6];
    char* histfile = new char[strlen(outFile)+6];
    char* statfile = new char[strlen(outFile)+6];
    strcpy(errorfile,outFile);
    strcat(errorfile,".errs");
    strcpy(histfile,outFile);
    strcat(histfile,".hist");
    strcpy(statfile,outFile);
    strcat(statfile,".stat");
    // Flush the old content of the files if they already exist.
    FILE* fp;
    fp = fopen(errorfile,"w");
    fclose(fp);
    fp = fopen(histfile,"w");
    fclose(fp);
    fp = fopen(statfile,"w");
    fprintf(fp,"Min\tMax\tVertexMean\tVertexSigma\tFaceMean\tFaceRMS\tHausdorff\tMSD\tMAD\tMedian\tDistVolume\t68percentile1\t68percentile2\t95percentile1\t95percentile2\t75percentile1\t75percentile2\n");
    fclose(fp);
    delete [] errorfile;
    delete [] histfile;
    delete [] statfile;
    
    if(batch)
    {
        int num;
        FILE* file = fopen(batchFile, "r");
        if(file == NULL)
        {
          fprintf(stdout,  "ERROR: couldn't open the batchFile %s\n",batchFile);
          return  -1;
        }
          
        if (fscanf(file, "%d", &num) != 1)
        {
          fprintf(stdout,  "ERROR reading the batchFile %s\n",batchFile);
          return  -1;
        }
        
        inFile1 = new char[500];
        inFile2 = new char[500];
        outFile = new char[500];
        
        for (int i=0; i < num; i++) 
        {
          if (fscanf(file, "%s", inFile1) != 1)
          {
                fprintf(stdout,  "ERROR reading the batchFile %s\n",batchFile);
                return  -1;
          }
          if (fscanf(file, "%s", inFile2) != 1)
          {
                fprintf(stdout,  "ERROR reading the batchFile %s\n",batchFile);
                return  -1;
          }
          if (fscanf(file, "%s", outFile) != 1)
          {
            fprintf(stdout,  "ERROR reading the batchFile %s\n",batchFile);
            return  -1;
          }
          fprintf(stdout,  "Comparing %s and %s ......\n",inFile1,inFile2);
          int res = compare(inFile1, inFile2, outFile, ss, freq, type, bins, absolute);
          if(res != 0)
            fprintf(stdout,  "failed.\n");
          else
            fprintf(stdout,  "completed.\n");
          
          delete [] inFile1;
          delete [] inFile2;
          delete [] outFile;
        }
    }
    else
    {
        fprintf(stdout,  "Comparing %s and %s ......\n",inFile1,inFile2);
        int res = compare(inFile1, inFile2, outFile, ss, freq, type, bins, absolute);
        if(res != 0)
          fprintf(stdout,  "failed.\n");
        else
                fprintf(stdout,  "completed.\n");
    }        
  }
  
  return 0;
}
