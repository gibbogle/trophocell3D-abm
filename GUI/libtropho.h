#ifndef LIBPARA32_H
#define LIBPARA32_H

#ifdef __cplusplus
extern "C" {
#endif
//
//
void execute(int *,char *, int *,char *, int *);
void simulate_step(int *);
void terminate_run(int *);
//void get_dimensions(int *,int *,int *);
void get_scene_dimensions(int *, int *, int *);
void get_scene(int *, int *);
void get_summary(int *);
//
//
#ifdef __cplusplus
}
#endif

#endif // LIBPARA32_H
